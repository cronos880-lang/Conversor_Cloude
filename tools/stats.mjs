import fs from 'fs'; import path from 'path';
import { parseSTL, buildMesh, components } from '../static/core/mesh.js';
const dir = process.argv[2];
for (const f of fs.readdirSync(dir).filter(f => /\.stl$/i.test(f)).sort()) {
  const buf = fs.readFileSync(path.join(dir, f));
  const t0 = performance.now();
  const m = buildMesh(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
  const c = components(m);
  console.log(`${f.padEnd(42)} tri=${String(m.nt).padStart(6)} v=${String(m.nv).padStart(6)} comp=${String(c.count).padStart(3)} open=${m.open} nm=${m.nonManifold} flip=${m.flipped} zero=${m.zeroArea} drop=${m.dropped} size=${m.max.map((v, i) => (v - m.min[i]).toFixed(1)).join('x')} vol=${m.volume.toFixed(0)} ${(performance.now() - t0).toFixed(0)}ms`);
}
