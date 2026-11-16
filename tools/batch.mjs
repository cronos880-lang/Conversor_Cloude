// Convierte todos los STL de una carpeta y valida cada STEP con Open CASCADE (python).
// node tools/batch.mjs carpeta_stl carpeta_salida [filtro]
import fs from 'fs'; import path from 'path'; import { execFileSync } from 'child_process';
import { parseSTL } from '../static/core/mesh.js';
import { convert } from '../static/core/convert.js';
const [, , dir, out, only] = process.argv; fs.mkdirSync(out, { recursive: true });
let okCount = 0, total = 0, faces = 0, tris = 0, secs = 0;
for (const f of fs.readdirSync(dir).filter(f => /\.stl$/i.test(f) && (!only || f.includes(only))).sort()) {
  total++; const name = f.replace(/\.stl$/i, ''); const buf = fs.readFileSync(path.join(dir, f));
  let line = name.slice(0, 30).padEnd(30);
  try {
    const t0 = performance.now();
    const r = convert(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)), { name });
    const s = (performance.now() - t0) / 1000; const file = path.join(out, name + '.step'); fs.writeFileSync(file, r.step);
    const rep = r.report; const ty = rep.types, cu = rep.curves;
    line += ` tri=${String(rep.triangles).padStart(6)} caras=${String(rep.faces).padStart(5)} [pl ${ty.plane || 0} ci ${ty.cylinder || 0} co ${ty.cone || 0} es ${ty.sphere || 0} to ${ty.torus || 0}] spl=${cu.bspline || 0} dev=${rep.edgeDeviation.toExponential(0)} ${s.toFixed(1)}s ${(r.step.length / 1024) | 0}KB`;
    let chk; try { chk = JSON.parse(execFileSync('python3', [new URL('./check_step.py', import.meta.url).pathname, file, String(rep.volume), String(rep.area)], { encoding: 'utf8', timeout: 120000 }).trim().split('\n').pop()); } catch (e) { chk = { error: 'lectura fallida' }; }
    const good = chk.valid && chk.solids >= 1 && Math.abs(chk['vol_err_%']) < 2 && Math.abs(chk['area_err_%']) < 2;
    line += ` | OCC ${good ? 'OK ' : 'MAL'} valid=${chk.valid} sol=${chk.solids} caras=${chk.faces} bad=${chk.bad_faces} vol=${chk['vol_err_%']}% area=${chk['area_err_%']}% ${chk.read_s}s ${chk.error || ''}`;
    if (good) okCount++; faces += rep.faces; tris += rep.triangles; secs += s;
  } catch (e) { line += ' ERROR ' + e.message; }
  console.log(line);
}
console.log(`TOTAL ${okCount}/${total} válidos · ${tris} triángulos -> ${faces} caras · ${secs.toFixed(0)}s`);
