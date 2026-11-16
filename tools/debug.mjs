import fs from 'fs';
import { parseSTL, buildMesh } from '../static/core/mesh.js';
import { segment } from '../static/core/segment.js';
const buf = fs.readFileSync(process.argv[2]);
const m = buildMesh(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
const { reg, prims } = segment(m);
// regiones planas pequeñas (candidatas a superficie curva no detectada)
const cnt = new Map(), cen = new Map();
for (let t = 0; t < m.nt; t++) { cnt.set(reg[t], (cnt.get(reg[t]) || 0) + 1); if (!cen.has(reg[t])) cen.set(reg[t], [0, 1, 2].map(j => (m.V[m.T[t * 3] * 3 + j] + m.V[m.T[t * 3 + 1] * 3 + j] + m.V[m.T[t * 3 + 2] * 3 + j]) / 3)); }
const small = [...cnt.entries()].filter(([r, c]) => prims[r].S.type === 'plane' && c <= 4);
console.log('planos pequeños:', small.length, 'ejemplos:', small.slice(0, 5).map(([r]) => cen.get(r).map(x => +x.toFixed(2)).join(',')).join(' | '));
if (process.argv[3]) { const c = process.argv[3].split(',').map(Number), rad = +process.argv[4] || 0.5; let n = 0; segment(m, { debug: p => Math.hypot(p[0] - c[0], p[1] - c[1], p[2] - c[2]) < rad && n++ < 6 }); }
