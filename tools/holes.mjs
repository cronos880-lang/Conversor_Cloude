import fs from 'fs'; import path from 'path';
import { parseSTL, buildMesh } from '../static/core/mesh.js';
const dir = process.argv[2];
for (const f of fs.readdirSync(dir).filter(f => /\.stl$/i.test(f)).sort()) {
  const buf = fs.readFileSync(path.join(dir, f));
  const m = buildMesh(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
  if (!m.open && !m.nonManifold) continue;
  const edges = new Map();
  for (let h = 0; h < m.nt * 3; h++) { const t = (h / 3) | 0, a = m.T[h], b = m.T[t * 3 + (h % 3 + 1) % 3]; const k = a < b ? a + ':' + b : b + ':' + a; (edges.get(k) || edges.set(k, []).get(k)).push([a, b, t]); }
  const P = i => [m.V[i * 3], m.V[i * 3 + 1], m.V[i * 3 + 2]].map(x => +x.toFixed(4));
  console.log(f, 'open', m.open, 'nm', m.nonManifold, 'diag', m.diag.toFixed(1));
  let n = 0;
  for (const [k, l] of edges) if (l.length !== 2 && n++ < 10) { const [a, b] = l[0]; console.log('   ', l.length, 'usos', P(a).join(','), '->', P(b).join(','), 'len', Math.hypot(...P(a).map((x, i) => x - P(b)[i])).toFixed(4), 'dirs', l.map(x => x[0] === a ? '+' : '-').join('')); }
}
