// node tools/one.mjs pieza.stl [salida.step]  -> conversión validada con traza de cada intento
import fs from 'fs'; import path from 'path';
import { parseSTL } from '../static/core/mesh.js';
import { convertValidated } from '../static/core/convert.js';
const init = (await import('../static/vendor/replicad_single.js')).default;
const wasmBinary = fs.readFileSync(new URL('../static/vendor/replicad_single.wasm', import.meta.url));
const oc = () => init({ wasmBinary, print: () => {}, printErr: () => {} });
const buf = fs.readFileSync(process.argv[2]); const t0 = performance.now();
try {
  const r = await convertValidated(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)), oc, { name: path.basename(process.argv[2]).replace(/\.stl$/i, ''), onStatus: s => process.env.QUIET || console.log('  ·', s), tolFactor: +(process.env.TOLF || 1), facetRegions: process.env.FACET ? (p => process.env.FACET.split(',').includes(p.S.type)) : process.env.FACETR ? ((p, i) => process.env.FACETR.split(',').map(Number).includes(i)) : process.env.LISTCYL ? ((p, i) => { if (p.S.type === 'cylinder') console.log('CIL', i, p.S.r.toFixed(3), p.S.a.map(x => x.toFixed(2)).join(','), p.S.p.map(x => x.toFixed(1)).join(',')); return false; }) : null, trace: c => console.log('    ', JSON.stringify(c)), maxRounds: +(process.env.ROUNDS || Infinity), dump: process.env.DUMP ? (r, step) => fs.writeFileSync(`${process.env.DUMP}_${r}.step`, step) : null });
  if (process.argv[3]) fs.writeFileSync(process.argv[3], r.step);
  console.log(JSON.stringify(r.report), ((performance.now() - t0) / 1000).toFixed(1) + 's');
} catch (e) { console.log('ERROR', typeof e, e && e.message || e); }
