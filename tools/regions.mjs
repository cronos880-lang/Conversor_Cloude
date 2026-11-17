// Regiones reconocidas (tras unificación y extrusiones): node tools/regions.mjs pieza.stl [n]
import fs from 'fs';
import { parseSTL, buildMesh, finishMesh } from '../static/core/mesh.js';
import { segment, mergeSameSurfaces } from '../static/core/segment.js';
import { findExtrusions } from '../static/core/extrude.js';
const b = fs.readFileSync(process.argv[2]);
let m = buildMesh(parseSTL(b.buffer.slice(b.byteOffset, b.byteOffset + b.length)));
if (m.volume < 0) { const T = m.T.slice(); for (let t = 0; t < m.nt; t++) { const x = T[t * 3 + 1]; T[t * 3 + 1] = T[t * 3 + 2]; T[t * 3 + 2] = x; } m = finishMesh(m.V, T, 0); }
const s = segment(m); const before = new Set(s.reg).size; mergeSameSurfaces(m, s); findExtrusions(m, s); mergeSameSurfaces(m, s);
const info = new Map();
for (let t = 0; t < m.nt; t++) { const r = s.reg[t]; const x = info.get(r) || { tri: 0, area: 0 }; x.tri++; x.area += m.A[t]; info.set(r, x); }
const rows = [...info].map(([r, x]) => ({ r, ...x, S: s.prims[r].S })).sort((a, b) => b.area - a.area);
const fmt = S => S.type === 'plane' ? '' : S.type === 'torus' ? `R ${S.R.toFixed(3)} r ${S.r.toFixed(3)}` : S.type === 'cone' ? `alfa ${(S.alpha * 57.3).toFixed(1)}` : S.r ? `r ${S.r.toFixed(3)}` : '';
console.log('regiones', before, '->', rows.length);
for (const x of rows.slice(0, +(process.argv[3] || 15))) console.log(x.S.type.padEnd(9), String(x.tri).padStart(6), 'tri', x.area.toFixed(2).padStart(10), 'mm2', fmt(x.S));
