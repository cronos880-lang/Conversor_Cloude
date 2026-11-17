"""Construye la versión descargable: python tools/build_entrega.py VERSION
Crea entrega/Conversor_Cloude_<VERSION>.zip con la aplicación web, el servidor local y los lanzadores."""
import json, sys, zipfile, hashlib, pathlib, datetime
root = pathlib.Path(__file__).resolve().parent.parent
version = sys.argv[1] if len(sys.argv) > 1 else 'C1'
name = f'Conversor_Cloude_{version}'
out = root / 'entrega' / f'{name}.zip'
files = [root / 'index.html'] + sorted(p for p in (root / 'static').rglob('*') if p.is_file())
manifest = {}
with zipfile.ZipFile(out, 'w', zipfile.ZIP_DEFLATED, compresslevel=9) as z:
    for p in files:
        data = p.read_bytes(); rel = 'app/' + p.relative_to(root).as_posix()
        z.writestr(f'{name}/{rel}', data); manifest[rel] = hashlib.sha256(data).hexdigest()
    for fn, crlf, bom in [('Iniciar_Conversor.bat', True, False), ('servidor_local.ps1', True, True), ('LEEME_PRIMERO.txt', True, False), ('iniciar_conversor.sh', False, False)]:
        text = (root / 'entrega' / fn).read_text(encoding='utf-8')
        if crlf: text = text.replace('\r\n', '\n').replace('\n', '\r\n')
        data = (b'\xef\xbb\xbf' if bom else b'') + text.encode('utf-8')
        info = zipfile.ZipInfo(f'{name}/{fn}', date_time=datetime.datetime.now().timetuple()[:6])
        info.compress_type = zipfile.ZIP_DEFLATED
        if fn.endswith('.sh'): info.external_attr = 0o755 << 16
        z.writestr(info, data)
    z.writestr(f'{name}/VERSION.json', json.dumps({'version': version, 'fecha': datetime.date.today().isoformat(), 'archivos': manifest}, indent=1, ensure_ascii=False))
print(out, out.stat().st_size, 'bytes')
