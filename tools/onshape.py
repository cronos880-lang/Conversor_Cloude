"""Verificación de STEP en Onshape real (API REST). Las claves las inyecta el proxy de la sesión.
python tools/onshape.py importar pieza.step [más.step ...]   -> importa en el documento de pruebas y resume
Guarda el id del documento de pruebas en tools/.onshape_doc.json (no contiene claves)."""
import json, sys, time, pathlib, requests
API = 'https://cad.onshape.com/api/v6'
H = {'Accept': 'application/json'}
STATE = pathlib.Path(__file__).with_name('.onshape_doc.json')

def call(method, path, **kw):
    for intento in range(5):
        r = requests.request(method, API + path, headers={**H, **kw.pop('headers', {})}, timeout=120, **kw)
        if r.status_code == 429:  # límite de llamadas: esperar y reintentar
            time.sleep(int(r.headers.get('Retry-After', '10'))); continue
        if r.status_code >= 400: raise RuntimeError(f'{method} {path}: HTTP {r.status_code} {r.text[:300]}')
        return r.json() if r.text else {}
    raise RuntimeError('Demasiados reintentos por límite de llamadas')

def documento():
    if STATE.exists():
        d = json.loads(STATE.read_text())
        try: call('GET', f"/documents/{d['did']}"); return d
        except RuntimeError: pass
    doc = call('POST', '/documents', json={'name': 'Pruebas conversor', 'isPublic': False})
    d = {'did': doc['id'], 'wid': doc['defaultWorkspace']['id']}
    STATE.write_text(json.dumps(d)); return d

def importar(path):
    d = documento(); p = pathlib.Path(path)
    with open(p, 'rb') as fh:
        tr = call('POST', f"/translations/d/{d['did']}/w/{d['wid']}", files={'file': (p.name, fh, 'application/octet-stream')},
                  data={'formatName': 'STEP', 'storeInDocument': 'false', 'flattenAssemblies': 'true', 'translate': 'true', 'createComposite': 'false', 'allowFaultyParts': 'true'})
    t0 = time.time()
    while tr.get('requestState') == 'ACTIVE' and time.time() - t0 < 900:
        time.sleep(3); tr = call('GET', f"/translations/{tr['id']}")
    return d, tr

if __name__ == '__main__':
    cmd, files = sys.argv[1], sys.argv[2:]
    if cmd == 'importar':
        for f in files:
            d, tr = importar(f)
            print(json.dumps({k: tr.get(k) for k in ('requestState', 'failureReason', 'resultElementIds', 'resultDocumentId', 'resultWorkspaceId')}, ensure_ascii=False))
