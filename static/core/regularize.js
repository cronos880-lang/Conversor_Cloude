// Regularización: relaciones exactas entre superficies vecinas (como en un diseño CAD).
// Pequeños errores de ajuste hacen que superficies que deberían ser paralelas, coaxiales o
// tangentes se crucen ligeramente cerca de su borde común, y Parasolid (Onshape) rechaza el
// sólido. Aquí se corrigen esos casos cuando la diferencia está dentro de la tolerancia.
import { dot, sub, add, mul, cross, norm, unit, dist } from './fit.js';

const RAD = Math.PI / 180;
const axisOf = S => S.a, originOf = S => S.p || S.apex || S.c;
const isRev = S => S.type === 'cylinder' || S.type === 'cone' || S.type === 'torus';

/** Distancia de un punto a la recta (o, a). */
const lineDist = (x, o, a) => { const v = sub(x, o); return norm(sub(v, mul(a, dot(v, a)))); };

export function regularize(mesh, seg) {
  const { V, T, adj, nt } = mesh, { reg, prims, tol } = seg;
  const R = prims.length, tris = prims.map(() => []), area = new Float64Array(R);
  for (let t = 0; t < nt; t++) { tris[reg[t]].push(t); area[reg[t]] += mesh.A[t]; }
  const live = i => tris[i].length && !prims[i].facet;
  const vertsOf = i => { const s = new Set(); for (const t of tris[i]) for (let j = 0; j < 3; j++) s.add(T[t * 3 + j]); return [...s].map(v => [V[v * 3], V[v * 3 + 1], V[v * 3 + 2]]); };
  // Solo se acepta un cambio si la región sigue sobre su superficie con algo de margen.
  const fitsAll = (i, S, k = 2) => vertsOf(i).every(p => Math.abs(dist(S, p)) <= k * tol);
  let changes = 0;
  const set = (i, S) => { if (fitsAll(i, S)) { prims[i] = { ...prims[i], S }; changes++; return true; } return false; };

  // Vecindad entre regiones.
  const nb = prims.map(() => new Set());
  for (let h = 0; h < nt * 3; h++) { if (adj[h] < 0) continue; const a = reg[(h / 3) | 0], b = reg[(adj[h] / 3) | 0]; if (a !== b) { nb[a].add(b); nb[b].add(a); } }

  // 1) Ejes de revolución casi paralelos a una dirección de referencia (normales de planos
  //    grandes y ejes de regiones mayores) -> exactamente paralelos.
  const order = prims.map((_, i) => i).filter(live).sort((a, b) => area[b] - area[a]);
  const refs = [];
  for (const i of order) {
    const S = prims[i].S;
    if (S.type === 'plane' && area[i] > 0.005 * mesh.area) refs.push(S.n);
    if (!isRev(S)) continue;
    let a = axisOf(S), best = null;
    for (const d of refs) { const c = Math.abs(dot(a, d)); if (c > Math.cos(1.5 * RAD) && c < 1 - 1e-12 && (!best || c > best.c)) best = { d, c }; }
    if (best) {
      const d = dot(a, best.d) < 0 ? mul(best.d, -1) : best.d, o = originOf(S);
      const S2 = S.type === 'cylinder' ? { ...S, a: d } : S.type === 'cone' ? { ...S, a: d } : { ...S, a: d };
      // Mantener el centro: el eje gira alrededor del punto medio de la región.
      const pts = vertsOf(i), m = pts.reduce((s, p) => add(s, mul(p, 1 / pts.length)), [0, 0, 0]);
      const h = dot(sub(m, o), a), pivot = add(o, mul(a, h));
      if (S.type === 'cylinder') S2.p = pivot;
      set(i, S2);
    }
    refs.push(prims[i].S.a);
  }

  // 2) Coaxialidad: superficies de revolución con ejes paralelos y casi coincidentes.
  const rev = order.filter(i => isRev(prims[i].S));
  for (let x = 0; x < rev.length; x++) for (let y = x + 1; y < rev.length; y++) {
    const i = rev[x], j = rev[y], Si = prims[i].S, Sj = prims[j].S;
    if (Math.abs(dot(Si.a, Sj.a)) < 1 - 1e-9) continue;
    const oi = originOf(Si), oj = originOf(Sj), off = lineDist(oj, oi, Si.a);
    if (off < 1e-9 || off > 2 * tol) continue;
    // Llevar el eje de la menor (j) al de la mayor (i).
    const v = sub(oj, oi), on = add(oi, mul(Si.a, dot(v, Si.a)));
    const S2 = { ...Sj }; if (Sj.type === 'cylinder') S2.p = on; else if (Sj.type === 'cone') S2.apex = add(on, sub(Sj.apex, oj)); else S2.c = on;
    if (Sj.type === 'cone') S2.apex = add(oi, mul(Si.a, dot(sub(Sj.apex, oi), Si.a)));
    set(j, S2);
  }

  // 3) Tangencias de redondeos (toros) con sus vecinos: cilindro coaxial o plano perpendicular al eje.
  for (const i of order) {
    const S = prims[i].S; if (S.type !== 'torus') continue;
    for (const j of nb[i]) {
      if (!live(j)) continue;
      const N = prims[j].S;
      if (N.type === 'cylinder' && Math.abs(dot(N.a, S.a)) > 1 - 1e-9 && lineDist(S.c, N.p, N.a) < 1e-6) {
        // Tubo tangente al cilindro: R - r = rc (por fuera) o R + r = rc (por dentro).
        for (const target of [N.r + S.r, N.r - S.r]) if (Math.abs(S.R - target) < 2 * tol && Math.abs(S.R - target) > 1e-12 && target > S.r) { set(i, { ...prims[i].S, R: target }); break; }
      } else if (N.type === 'plane' && Math.abs(dot(N.n, S.a)) > 1 - 1e-9) {
        // Tubo tangente al plano: la altura del centro del tubo dista r del plano.
        const cur = prims[i].S, hc = dot(N.n, cur.c) - N.d;
        for (const target of [cur.r, -cur.r]) if (Math.abs(hc - target) < 2 * tol && Math.abs(hc - target) > 1e-12) { set(i, { ...cur, c: add(cur.c, mul(N.n, target - hc)) }); break; }
      }
    }
  }

  // 4) Cilindros de redondeo tangentes a planos paralelos a su eje: distancia eje-plano = r.
  for (const i of order) {
    const S = prims[i].S; if (S.type !== 'cylinder') continue;
    for (const j of nb[i]) {
      if (!live(j) || prims[j].S.type !== 'plane') continue;
      const P = prims[j].S; if (Math.abs(dot(P.n, prims[i].S.a)) > 1e-9) continue;
      const cur = prims[i].S, d0 = dot(P.n, cur.p) - P.d;
      for (const target of [cur.r, -cur.r]) if (Math.abs(d0 - target) < 2 * tol && Math.abs(d0 - target) > 1e-12) { set(i, { ...cur, p: add(cur.p, mul(P.n, target - d0)) }); break; }
    }
  }
  return changes;
}
void cross; void unit;
