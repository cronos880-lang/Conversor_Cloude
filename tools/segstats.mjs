import fs from 'fs'; import path from 'path';
import { parseSTL, buildMesh } from '../static/core/mesh.js';
import { segment } from '../static/core/segment.js';
const dir = process.argv[2], only = process.argv[3];
for (const f of fs.readdirSync(dir).filter(f => /\.stl$/i.test(f) && (!only || f.includes(only))).sort()) {
  const buf = fs.readFileSync(path.join(dir, f));
  const m = buildMesh(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
  const t0 = performance.now();
  const { reg, prims, tol } = segment(m);
  const dt = performance.now() - t0;
  const count = {}, area = {};
  for (const p of prims) count[p.S.type] = (count[p.S.type] || 0) + 1;
  for (let t = 0; t < m.nt; t++) { const ty = prims[reg[t]].S.type; area[ty] = (area[ty] || 0) + m.A[t]; }
  const ord = ['plane', 'cylinder', 'cone', 'sphere', 'torus'];
  console.log(`${f.slice(0, 34).padEnd(34)} tri=${String(m.nt).padStart(6)} reg=${String(prims.length).padStart(5)} ${ord.map(k => `${k.slice(0, 3)}=${count[k] || 0}(${(100 * (area[k] || 0) / m.area).toFixed(0)}%)`).join(' ')} tol=${tol.toExponential(1)} ${(dt / 1000).toFixed(1)}s`);
}
