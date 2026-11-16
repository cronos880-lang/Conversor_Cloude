import fs from 'fs'; import path from 'path';
import { parseSTL, buildMesh } from '../static/core/mesh.js';
import { segment } from '../static/core/segment.js';
const dir = process.argv[2];
for (const f of fs.readdirSync(dir).filter(f => /\.stl$/i.test(f)).sort()) {
  const buf = fs.readFileSync(path.join(dir, f)); const head = buf.slice(0, 5).toString();
  const m = buildMesh(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
  const { reg, prims, tol } = segment(m); const area = prims.map(() => 0); for (let t = 0; t < m.nt; t++) area[reg[t]] += m.A[t];
  const q = (type) => { const xs = prims.map((p, i) => [p.res, area[i], p.S.type]).filter(x => (type === 'plane') === (x[2] === 'plane') && x[1] > 0).sort((a, b) => a[0] - b[0]); const tot = xs.reduce((s, x) => s + x[1], 0); let acc = 0; const out = []; for (const p of [0.5, 0.9, 1]) { acc = 0; for (const x of xs) { acc += x[1]; if (acc >= p * tot - 1e-12) { out.push(x[0].toExponential(0)); break; } } } return out.join('/'); };
  console.log(f.slice(0, 32).padEnd(32), head === 'solid' ? 'ASCII?' : 'bin   ', 'tol', tol.toFixed(3), 'planos p50/p90/max', q('plane').padEnd(20), 'curvas', q('curved'));
}
