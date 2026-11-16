// node tools/facedump.mjs pieza.stl x,y,z  -> describe la cara cuyo centro de área está más cerca del punto
import fs from 'fs';
import { parseSTL } from '../static/core/mesh.js';
import { convert } from '../static/core/convert.js';
const buf = fs.readFileSync(process.argv[2]); const c = process.argv[3].split(',').map(Number); const type = process.argv[4];
const r = convert(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
const { mesh: m, seg, brep } = r; const cen = brep.faces.map(() => [0, 0, 0, 0]);
for (let t = 0; t < m.nt; t++) { const f = cen[seg.reg[t]]; for (let j = 0; j < 3; j++) f[j] += m.A[t] * (m.V[m.T[t * 3] * 3 + j] + m.V[m.T[t * 3 + 1] * 3 + j] + m.V[m.T[t * 3 + 2] * 3 + j]) / 3; f[3] += m.A[t]; }
let best = -1, bd = Infinity;
cen.forEach((f, i) => { if (!f[3] || (type && brep.faces[i].S.type !== type)) return; const d = Math.hypot(f[0] / f[3] - c[0], f[1] / f[3] - c[1], f[2] / f[3] - c[2]); if (d < bd) { bd = d; best = i; } });
const f = brep.faces[best]; const fx = v => v.map(x => +x.toFixed(4)).join(',');
console.log('cara', best, 'dist', bd.toFixed(3), 'area', cen[best][3].toFixed(4), JSON.stringify(f.S, (k, v) => typeof v === 'number' ? +v.toFixed(5) : v), 'sign', f.sign, 'res', seg.prims[best].res.toExponential(1));
f.loops.forEach((loop, i) => { console.log(' bucle', i); for (const ref of loop) { const e = brep.edges[ref.edge]; const cu = e.curve; console.log('   ', ref.forward ? '+' : '-', 'e' + ref.edge, cu.type, 'n=' + e.n, 'dev=' + e.dev.toExponential(1), fx(brep.vertices[e.v0]), '->', fx(brep.vertices[e.v1]), cu.type === 'circle' ? `r=${cu.r.toFixed(4)} c=${fx(cu.c)} a=${fx(cu.a)} sweep=${(cu.sweep * 57.3).toFixed(1)}` : ''); } });
