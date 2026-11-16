// Pruebas del núcleo con el motor WebAssembly real: node tests/core_test.mjs
import fs from 'fs';
import { convertValidated } from '../static/core/convert.js';
const init = (await import('../static/vendor/replicad_single.js')).default;
const newOc = () => init({ wasmBinary: fs.readFileSync(new URL('../static/vendor/replicad_single.wasm', import.meta.url)), print: () => {}, printErr: () => {} });

const flat = tris => new Float32Array(tris.flat(2));
function tube(n, R, r, h) {
  const T = [], P = (rad, i, z) => [rad * Math.cos(2 * Math.PI * i / n), rad * Math.sin(2 * Math.PI * i / n), z];
  for (let i = 0; i < n; i++) {
    const j = (i + 1) % n;
    T.push([P(R, i, 0), P(R, j, 0), P(R, j, h)], [P(R, i, 0), P(R, j, h), P(R, i, h)], [P(r, i, 0), P(r, j, h), P(r, j, 0)], [P(r, i, 0), P(r, i, h), P(r, j, h)]);
    T.push([P(r, i, h), P(R, i, h), P(R, j, h)], [P(r, i, h), P(R, j, h), P(r, j, h)], [P(r, i, 0), P(R, j, 0), P(R, i, 0)], [P(r, i, 0), P(r, j, 0), P(R, j, 0)]);
  }
  return T;
}
function sphere(nu, nv, R) {
  const T = [], P = (i, j) => { const th = Math.PI * j / nv, ph = 2 * Math.PI * (i % nu) / nu; return j === 0 ? [0, 0, R] : j === nv ? [0, 0, -R] : [R * Math.sin(th) * Math.cos(ph), R * Math.sin(th) * Math.sin(ph), R * Math.cos(th)]; };
  for (let i = 0; i < nu; i++) for (let j = 0; j < nv; j++) { const a = P(i, j), b = P(i + 1, j), c = P(i + 1, j + 1), d = P(i, j + 1); if (j > 0) T.push([a, d, b]); if (j < nv - 1) T.push([b, d, c]); }
  return T;
}
function box(o, s, flip) {
  const v = []; for (let i = 0; i < 8; i++) v.push([o[0] + (i & 1 ? s : 0), o[1] + (i & 2 ? s : 0), o[2] + (i & 4 ? s : 0)]);
  const T = [];
  for (const [a, b, c, d] of [[0, 2, 3, 1], [4, 5, 7, 6], [0, 1, 5, 4], [2, 6, 7, 3], [0, 4, 6, 2], [1, 3, 7, 5]]) for (const t of [[a, b, c], [a, c, d]]) T.push((flip ? [t[0], t[2], t[1]] : t).map(i => v[i]));
  return T;
}
const openBox = box([0, 0, 0], 5, false).slice(0, 11);

const cases = [
  ['caja', box([0, 0, 0], 10, false), { faces: 6, types: { plane: 6 }, volume: 1000 }],
  ['tubo de 128 lados', tube(128, 20, 10, 30), { faces: 4, types: { cylinder: 2, plane: 2 }, volume: Math.PI * 300 * 30 }],
  ['tubo en pulgadas', tube(64, 2, 1, 3), { faces: 4, types: { cylinder: 2, plane: 2 }, volume: Math.PI * 3 * 3 * 25.4 ** 3, scale: 25.4 }],
  ['esfera', sphere(90, 45, 10), { faces: 2, types: { sphere: 2 }, volume: 4 / 3 * Math.PI * 1000 }],
  ['caja con cavidad', [...box([0, 0, 0], 10, false), ...box([3, 3, 3], 4, true)], { faces: 12, types: { plane: 12 }, volume: 936 }],
  ['dos piezas', [...box([0, 0, 0], 10, false), ...box([20, 0, 0], 5, false)], { faces: 12, types: { plane: 12 }, volume: 1125 }],
  ['caja con un hueco (reparación)', openBox, { faces: 6, types: { plane: 6 }, volume: 125 }],
];
let failed = 0;
for (const [name, tris, want] of cases) {
  try {
    const r = convertValidated(flat(tris), await newOc(), { scale: want.scale || 1, name });
    const rep = r.report, errs = [];
    if (rep.faces !== want.faces) errs.push(`caras ${rep.faces} ≠ ${want.faces}`);
    for (const k of Object.keys(want.types)) if (rep.types[k] !== want.types[k]) errs.push(`${k} ${rep.types[k]} ≠ ${want.types[k]}`);
    if (Math.abs(rep.checkedVolume - want.volume) > 0.005 * want.volume) errs.push(`volumen ${rep.checkedVolume} ≠ ${want.volume}`);
    if (!r.step.startsWith('ISO-10303-21;')) errs.push('cabecera STEP');
    console.log(errs.length ? 'FALLA' : 'PASA ', name, errs.join('; '));
    if (errs.length) failed++;
  } catch (e) { console.log('FALLA', name, e.message || e); failed++; }
}
try { convertValidated(new Float32Array([0, 0, 0, 1, 0, 0, 0, 1, 0]), await newOc(), {}); console.log('FALLA triángulo suelto: debía rechazarse'); failed++; } catch { console.log('PASA  rechazo de malla abierta no reparable'); }
console.log(failed ? `${failed} pruebas fallidas` : 'Todas las pruebas superadas');
process.exit(failed ? 1 : 0);
