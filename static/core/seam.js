// Uniones exactas entre facetas y superficies reconocidas.
// Las facetas (zonas sin reconocer) siguen los vértices de la malla, que quedan hasta una
// tolerancia fuera de la superficie reconocida vecina: Parasolid (Onshape) ve entonces bordes
// separados de sus caras y llega a descartar caras enteras. Aquí, los triángulos de faceta que
// tocan una superficie reconocida pasan a ser caras sueltas, sus vértices compartidos se llevan
// exactamente sobre esa superficie y el plano de cada triángulo se recalcula (un triángulo
// sigue siendo plano aunque se muevan sus vértices).
import { snapPoint } from './curves.js';
import { cross, sub, norm, dot, mul } from './fit.js';

export function exactFacetJoins(mesh, seg, maxMoveFactor = 3) {
  const { V, T, nt, nv } = mesh, { reg, prims, tol } = seg;
  const P = i => [V[i * 3], V[i * 3 + 1], V[i * 3 + 2]];
  // Regiones que tocan cada vértice.
  const touch = Array.from({ length: nv }, () => null);
  for (let t = 0; t < nt; t++) for (let k = 0; k < 3; k++) { const v = T[t * 3 + k]; (touch[v] || (touch[v] = new Set())).add(reg[t]); }
  const mixed = v => { const s = touch[v]; if (!s) return false; let f = false, a = false; for (const r of s) { if (prims[r].facet) f = true; else a = true; } return f && a; };
  // 1) Triángulos de faceta con algún vértice en una unión: región propia.
  const count = new Map(); for (let t = 0; t < nt; t++) count.set(reg[t], (count.get(reg[t]) || 0) + 1);
  for (let t = 0; t < nt; t++) {
    const r = reg[t]; if (!prims[r].facet || count.get(r) === 1) continue;
    if (![0, 1, 2].some(k => mixed(T[t * 3 + k]))) continue;
    count.set(r, count.get(r) - 1);
    const id = prims.length; prims.push({ ...prims[r], facet: true, isolated: true }); reg[t] = id; count.set(id, 1);
  }
  // Recalcular contactos tras separar.
  for (let v = 0; v < nv; v++) touch[v] = null;
  for (let t = 0; t < nt; t++) for (let k = 0; k < 3; k++) { const v = T[t * 3 + k]; (touch[v] || (touch[v] = new Set())).add(reg[t]); }
  // 2) Vértices de la unión sobre las superficies reconocidas que tocan.
  const moved = new Map();
  for (let v = 0; v < nv; v++) {
    if (!mixed(v)) continue;
    const surf = [...touch[v]].filter(r => !prims[r].facet).map(r => prims[r].S);
    const p = P(v), q = snapPoint(p, surf, maxMoveFactor * tol);
    if (q !== p && norm(sub(q, p)) > 1e-12) moved.set(v, q);
  }
  // 3) Plano de cada triángulo de faceta aislado a partir de sus vértices (movidos o no).
  const at = v => moved.get(v) || P(v);
  for (let t = 0; t < nt; t++) {
    const r = reg[t]; if (!prims[r].isolated) continue;
    const [a, b, c] = [0, 1, 2].map(k => at(T[t * 3 + k])), n0 = cross(sub(b, a), sub(c, a)), l = norm(n0);
    if (!(l > 1e-14)) continue;
    let n = mul(n0, 1 / l); if (dot(n, prims[r].S.n) < 0) n = mul(n, -1);
    prims[r] = { ...prims[r], S: { type: 'plane', n, d: dot(n, a) } };
  }
  return moved;
}
