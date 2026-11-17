// Desviación de aristas por pareja de superficies: node tools/edgedev.mjs pieza.stl
import fs from 'fs';
import { parseSTL } from '../static/core/mesh.js';
import { convert } from '../static/core/convert.js';
const buf = fs.readFileSync(process.argv[2]);
const r = convert(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
const st = {};
for (const c of r.brep.chainInfo) { const k = `${c.A}/${c.B} ${c.kind}`; (st[k] = st[k] || []).push(c); }
const edgeDev = new Map(); r.brep.edges.forEach(e => {});
for (const [k, list] of Object.entries(st)) console.log(k.padEnd(32), list.length);
const bad = r.brep.edges.filter(e => !e.poly && e.dev > 0.1).sort((a, b) => b.dev - a.dev).slice(0, 6);
for (const e of bad) console.log('dev', e.dev.toFixed(3), e.curve.type, 'n', e.n, r.brep.vertices[e.v0].map(x => +x.toFixed(2)).join(','), '->', r.brep.vertices[e.v1].map(x => +x.toFixed(2)).join(','));
import { dist } from '../static/core/fit.js';
for (const e of bad) {
  const p = r.brep.vertices[e.v0], q = r.brep.vertices[e.v1], m = p.map((x, i) => (x + q[i]) / 2);
  const SA = r.seg.prims[e.A].S, SB = r.seg.prims[e.B].S;
  console.log(SA.type, SB.type, 'dA', [p, m, q].map(x => dist(SA, x).toFixed(4)).join('/'), 'dB', [p, m, q].map(x => dist(SB, x).toFixed(4)).join('/'), Array.isArray(SA.d) ? 'dirA ' + SA.d.map(x => x.toFixed(3)) : '', Array.isArray(SB.d) ? 'dirB ' + SB.d.map(x => x.toFixed(3)) : '');
}
{
  const e = bad[0]; const A = r.seg.prims[e.A].S.type === 'plane' ? e.A : e.B; const ts = []; for (let t = 0; t < r.mesh.nt; t++) if (r.seg.reg[t] === A) ts.push(t);
  const P = i => [0, 1, 2].map(k => r.mesh.V[i * 3 + k]);
  console.log('plano', A, JSON.stringify(r.seg.prims[A]).slice(0, 200), 'tris', ts.length, 'dist vértices', ts.flatMap(t => [0, 1, 2].map(j => dist(r.seg.prims[A].S, P(r.mesh.T[t * 3 + j])).toFixed(3))).join(' '));
}
