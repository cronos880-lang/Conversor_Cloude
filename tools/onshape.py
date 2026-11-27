"""Verificación de STEP en Onshape real (API REST). Las claves las inyecta el proxy de la sesión.
python tools/onshape.py verificar pieza.step [más.step ...]  -> importa, mide y borra lo importado
python tools/onshape.py limpiar                             -> borra del documento de pruebas todo lo importado
El documento de pruebas ("Pruebas conversor", público por ser cuenta gratuita) se guarda en tools/.onshape_doc.json."""
import json, sys, time, pathlib, requests
API = 'https://cad.onshape.com/api/v6'
STATE = pathlib.Path(__file__).with_name('.onshape_doc.json')
BASE_ELEMENTS = {'Part Studio 1', 'Assembly 1', 'BOM : Assembly 1'}

def call(method, path, **kw):
    for intento in range(6):
        r = requests.request(method, API + path, headers={'Accept': 'application/json'}, timeout=300, **kw)
        if r.status_code == 429:  # límite de llamadas: esperar y reintentar
            time.sleep(int(r.headers.get('Retry-After', '15'))); continue
        if r.status_code >= 400: raise RuntimeError(f'{method} {path}: HTTP {r.status_code} {r.text[:300]}')
        return r.json() if r.text else {}
    raise RuntimeError('Demasiados reintentos por límite de llamadas')

def documento():
    if STATE.exists():
        d = json.loads(STATE.read_text())
        try: call('GET', f"/documents/{d['did']}"); return d
        except RuntimeError: pass
    doc = call('POST', '/documents', json={'name': 'Pruebas conversor', 'isPublic': True})
    d = {'did': doc['id'], 'wid': doc['defaultWorkspace']['id']}
    STATE.write_text(json.dumps(d)); return d

def fs_value(v):
    """Convierte un valor FeatureScript serializado en Python."""
    t = v.get('btType', '')
    if 'Map' in t and isinstance(v.get('value'), list):
        return {fs_value(e['key']): fs_value(e['value']) for e in v['value']}
    if 'Array' in t: return [fs_value(x) for x in v['value']]
    return v.get('value')

SCRIPT = '''function(context is Context, queries) {
  var st = getFeatureStatus(context, makeId("%s"));
  var solids = qBodyType(qEverything(EntityType.BODY), BodyType.SOLID);
  var tipos = {};
  var caras = [];
  for (var f in evaluateQuery(context, qOwnedByBody(solids, EntityType.FACE))) {
    var k = try silent(toString(evSurfaceDefinition(context, { "face" : f }).surfaceType)); if (k == undefined) k = "ERROR";
    tipos[k] = (tipos[k] == undefined ? 0 : tipos[k]) + 1;
    if (%s) {
      var c = try silent(evApproximateCentroid(context, { "entities" : f })); var ar = try silent(evArea(context, { "entities" : f }) / (millimeter ^ 2));
      caras = append(caras, [k, ar, c == undefined ? undefined : c[0] / millimeter, c == undefined ? undefined : c[1] / millimeter, c == undefined ? undefined : c[2] / millimeter]);
    }
  }
  var vol = 0;
  for (var b in evaluateQuery(context, solids)) { var vb = try silent(evVolume(context, { "entities" : b }) / (millimeter ^ 3)); if (vb != undefined) vol += vb; }
  return { "estado" : toString(st.statusType), "motivo" : toString(st.statusEnum), "solidos" : size(evaluateQuery(context, solids)), "volumen" : vol, "caras_detalle" : caras,
           "superficies" : size(evaluateQuery(context, qBodyType(qEverything(EntityType.BODY), BodyType.SHEET))), "tipos" : tipos };
}'''

def verificar(path, borrar=False, detalle=False):
    d = documento(); p = pathlib.Path(path); t0 = time.time()
    with open(p, 'rb') as fh:
        blob = call('POST', f"/blobelements/d/{d['did']}/w/{d['wid']}", files={'file': (p.stem + '.step', fh, 'application/step')},
                    data={'encodedFilename': p.stem + '.step', 'translate': 'true', 'formatName': 'STEP', 'allowFaultyParts': 'true', 'flattenAssemblies': 'true', 'unit': 'millimeter'})
    # La cuenta gratuita tiene un cupo anual de llamadas: esperar según el tamaño y consultar pocas veces.
    espera = min(240, 8 + p.stat().st_size / 150000); tr = {'requestState': 'ACTIVE'}
    for intento in range(8):
        time.sleep(espera); tr = call('GET', f"/translations/{blob['translationId']}")
        if tr.get('requestState') != 'ACTIVE': break
        espera = min(300, espera * 1.6)
    out = {'pieza': p.stem, 'traduccion': tr.get('requestState'), 'segundos': round(time.time() - t0)}
    creados = [blob['id']] + list(tr.get('resultElementIds') or [])
    try:
        if tr.get('requestState') != 'DONE': out['fallo'] = tr.get('failureReason'); return out
        ps = (tr.get('resultElementIds') or [None])[0]
        base = f"/partstudios/d/{d['did']}/w/{d['wid']}/e/{ps}"
        feats = call('GET', base + '/features')['features']
        imp = [f['featureId'] for f in feats if f.get('featureType') == 'importForeign']
        res = fs_value(call('POST', base + '/featurescript', json={'script': SCRIPT % (imp[0], 'true' if detalle else 'false'), 'queries': {}})['result'])
        out.update(res); out['caras'] = sum((res.get('tipos') or {}).values())
        return out
    finally:
        if borrar:
            for e in creados:
                try: call('DELETE', f"/elements/d/{d['did']}/w/{d['wid']}/e/{e}")
                except RuntimeError: pass

def limpiar():
    d = documento()
    for e in call('GET', f"/documents/d/{d['did']}/w/{d['wid']}/elements"):
        if e['name'] not in BASE_ELEMENTS: call('DELETE', f"/elements/d/{d['did']}/w/{d['wid']}/e/{e['id']}")
    ps = [e for e in call('GET', f"/documents/d/{d['did']}/w/{d['wid']}/elements") if e['name'] == 'Part Studio 1']
    for e in ps:
        base = f"/partstudios/d/{d['did']}/w/{d['wid']}/e/{e['id']}"
        for f in call('GET', base + '/features').get('features', []): call('DELETE', base + '/features/featureid/' + f['featureId'])

if __name__ == '__main__':
    cmd, files = sys.argv[1], sys.argv[2:]
    if cmd == 'limpiar': limpiar(); print('documento limpio')
    elif cmd in ('verificar', 'detalle'):
        for f in files:
            r = verificar(f, detalle=cmd == 'detalle')
            if cmd == 'detalle':
                pathlib.Path(f).with_suffix('.onshape.json').write_text(json.dumps(r.pop('caras_detalle', []))); r['detalle'] = str(pathlib.Path(f).with_suffix('.onshape.json'))
            else: r.pop('caras_detalle', None)
            print(json.dumps(r, ensure_ascii=False), flush=True)
