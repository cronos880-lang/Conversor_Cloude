// Regularización: relaciones exactas entre superficies vecinas (como en un diseño CAD).
// Pequeños errores de ajuste hacen que superficies que deberían ser paralelas, coaxiales o
// tangentes se crucen ligeramente cerca de su borde común, y Parasolid (Onshape) rechaza el
// sólido. Aquí se corrigen esos casos cuando la diferencia está dentro de la tolerancia.
import { dot, sub, add, mul, cross, norm, unit, dist, perp, solve, cylinderWithAxis, normalAt, fitPlane } from './fit.js';

const RAD = Math.PI / 180;
const axisOf = S => S.a, originOf = S => S.p || S.apex || S.c;
const isRev = S => S.type === 'cylinder' || S.type === 'cone' || S.type === 'torus';

/** Distancia de un punto a la recta (o, a). */
const lineDist = (x, o, a) => { const v = sub(x, o); return norm(sub(v, mul(a, dot(v, a)))); };


/**
 * Cilindro o cono con la dirección del eje fijada (p. ej. un taladro con desmoldeo): ajusta
 * centro, radio y conicidad por mínimos cuadrados (Gauss-Newton).
 */
function revolutionWithAxis(P, d) {
  const C0 = cylinderWithAxis(P, d); if (!C0) return null;
  const e1 = perp(d), e2 = cross(d, e1), o = C0.p;
  const loc = P.map(p => { const v = sub(p, o); return [dot(v, e1), dot(v, e2), dot(v, d)]; });
  let x = [0, 0, C0.r, 0];
  for (let it = 0; it < 30; it++) {
    const M = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], g = [0, 0, 0, 0];
    for (const [u, v, h] of loc) {
      const du = u - x[0], dv = v - x[1], rho = Math.hypot(du, dv) || 1e-12, res = rho - (x[2] + x[3] * h);
      const J = [-du / rho, -dv / rho, -1, -h];
      for (let i = 0; i < 4; i++) { g[i] -= J[i] * res; for (let j = 0; j < 4; j++) M[i][j] += J[i] * J[j]; }
    }
    for (let i = 0; i < 4; i++) M[i][i] += 1e-12;
    const st = solve(M, g); if (!st) return null;
    x = x.map((xi, i) => xi + st[i]);
    if (Math.hypot(...st) < 1e-12) break;
  }
  const c = add(o, add(mul(e1, x[0]), mul(e2, x[1])));
  if (Math.abs(x[3]) < 1e-4) return { type: 'cylinder', p: c, a: d, r: x[2] };
  // Cono: el vértice está donde el radio se anula; el eje apunta hacia donde se ensancha.
  const hApex = -x[2] / x[3], apex = add(c, mul(d, hApex));
  return { type: 'cone', apex, a: x[3] > 0 ? d : mul(d, -1), alpha: Math.atan(Math.abs(x[3])) };
}

export function regularize(mesh, seg) {
  const { V, T, adj, nt } = mesh, { reg, prims, tol } = seg;
  const R = prims.length, tris = prims.map(() => []), area = new Float64Array(R);
  for (let t = 0; t < nt; t++) { tris[reg[t]].push(t); area[reg[t]] += mesh.A[t]; }
  const live = i => tris[i].length && !prims[i].facet;
  const vertsOf = i => { const s = new Set(); for (const t of tris[i]) for (let j = 0; j < 3; j++) s.add(T[t * 3 + j]); return [...s].map(v => [V[v * 3], V[v * 3 + 1], V[v * 3 + 2]]); };
  // Solo se acepta un cambio si la región sigue sobre su superficie con algo de margen.
  // Si el ajuste original ya se separaba más de 2·tol, se admite un empeoramiento pequeño.
  const worst = (i, S) => vertsOf(i).reduce((m, p) => Math.max(m, Math.abs(dist(S, p))), 0);
  const fitsAll = (i, S, k = 2) => worst(i, S) <= Math.max(k * tol, Math.min(4 * tol, worst(i, prims[i].S) + 0.25 * tol));
  let changes = 0;
  const set = (i, S) => { if (fitsAll(i, S)) { prims[i] = { ...prims[i], S }; changes++; return true; } return false; };

  // Vecindad entre regiones.
  const nb = prims.map(() => new Set());
  for (let h = 0; h < nt * 3; h++) { if (adj[h] < 0) continue; const a = reg[(h / 3) | 0], b = reg[(adj[h] / 3) | 0]; if (a !== b) { nb[a].add(b); nb[b].add(a); } }

  // 0) Planos casi paralelos o perpendiculares a planos mayores (< 0,5°): exactamente.
  if (!globalThis.process?.env?.NOPLANE) {
    const big = [];
    for (const i of prims.map((_, k) => k).filter(live).sort((a, b) => area[b] - area[a])) {
      const S = prims[i].S; if (S.type !== 'plane') continue;
      let n = S.n;
      for (const d of big) {
        const c = dot(n, d);
        if (Math.abs(c) > Math.cos(0.5 * RAD) && Math.abs(c) < 1 - 1e-15) n = c < 0 ? mul(d, -1) : d;
        else if (Math.abs(c) < Math.sin(0.5 * RAD) && Math.abs(c) > 1e-15) n = unit(sub(n, mul(d, c)));
      }
      if (n !== S.n) {
        const pts = vertsOf(i), m = pts.reduce((s2, p) => add(s2, mul(p, 1 / pts.length)), [0, 0, 0]);
        set(i, { ...S, n, d: dot(n, m) });
      }
      if (area[i] > 0.005 * mesh.area && big.length < 12 && !big.some(d => Math.abs(dot(d, prims[i].S.n)) > 1 - 1e-12)) big.push(prims[i].S.n);
    }
  }

  // 1) Ejes de revolución casi paralelos a una dirección de referencia (normales de planos
  //    grandes y ejes de regiones mayores) -> exactamente paralelos.
  const order = prims.map((_, i) => i).filter(live).sort((a, b) => area[b] - area[a]);
  const refs = [], planeRefs = [];
  for (const i of order) {
    const S = prims[i].S;
    if (S.type === 'plane' && area[i] > 0.005 * mesh.area) { refs.push(S.n); planeRefs.push(S.n); }
    if (!isRev(S)) continue;
    let a = axisOf(S), best = null;
    for (const d of refs) { const c = Math.abs(dot(a, d)); if (c > Math.cos(1.5 * RAD) && c < 1 - 1e-12 && (!best || c > best.c)) best = { d, c, plane: planeRefs.includes(d) }; }
    if (S.type === 'cylinder' && !(best && best.plane) && !globalThis.process?.env?.NOCONE) {
      // Eje inclinado unos grados respecto a una dirección principal: suele ser un taladro o
      // saliente con ángulo de desmoldeo (cono) que se ajustó como cilindro torcido.
      let near = null;
      for (const d of planeRefs) { const c = Math.abs(dot(a, d)); if (c > Math.cos(6 * RAD) && c <= Math.cos(1.5 * RAD) && (!near || c > near.c)) near = { d, c }; }
      if (near) { const S2 = revolutionWithAxis(vertsOf(i), near.d); if (S2 && set(i, S2)) { refs.push(near.d); continue; } }
    }
    if (best) {
      const d = dot(a, best.d) < 0 ? mul(best.d, -1) : best.d, o = originOf(S);
      const S2 = S.type === 'cylinder' ? { ...S, a: d } : S.type === 'cone' ? { ...S, a: d } : { ...S, a: d };
      // Mantener el centro: el eje gira alrededor del punto medio de la región.
      const pts = vertsOf(i), m = pts.reduce((s, p) => add(s, mul(p, 1 / pts.length)), [0, 0, 0]);
      const h = dot(sub(m, o), a), pivot = add(o, mul(a, h));
      if (S.type === 'cylinder') S2.p = pivot;
      set(i, S2);
    }
    // Eje casi contenido en la dirección de un plano grande (perpendicular a su normal): exacto.
    if (prims[i].S.type === 'cylinder') for (const d of planeRefs) {
      const C = prims[i].S, c = dot(C.a, d);
      if (Math.abs(c) < 1e-12 || Math.abs(c) > Math.sin(1.5 * RAD)) continue;
      const pts = vertsOf(i), m = pts.reduce((s2, p) => add(s2, mul(p, 1 / pts.length)), [0, 0, 0]);
      const pivot = add(C.p, mul(C.a, dot(sub(m, C.p), C.a)));
      set(i, { ...C, a: unit(sub(C.a, mul(d, c))), p: pivot });
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
  //    El plano puede estar girado unas décimas de grado: el eje del toro pasa a ser su normal.
  for (const i of order) {
    if (prims[i].S.type !== 'torus') continue;
    for (const j of nb[i]) {
      if (!live(j)) continue;
      const N = prims[j].S, S = prims[i].S;
      if (N.type === 'cylinder' && Math.abs(dot(N.a, S.a)) > Math.cos(1.5 * RAD) && lineDist(S.c, N.p, N.a) < 3 * tol) {
        // Tubo tangente al cilindro: R - r = rc (por fuera) o R + r = rc (por dentro), con el
        // mismo eje que el cilindro.
        const a = dot(N.a, S.a) < 0 ? mul(N.a, -1) : N.a, c = add(N.p, mul(N.a, dot(sub(S.c, N.p), N.a)));
        for (const target of [N.r + S.r, N.r - S.r]) if (Math.abs(S.R - target) < 3 * tol && target > S.r) {
          if (Math.abs(S.R - target) > 1e-12 || norm(sub(c, S.c)) > 1e-12 || Math.abs(dot(a, S.a)) < 1 - 1e-12) set(i, { ...S, a, c, R: target });
          break;
        }
      } else if (N.type === 'plane' && Math.abs(dot(N.n, S.a)) > Math.cos(1.5 * RAD)) {
        // Tubo tangente al plano: la altura del centro del tubo dista r del plano.
        const a = dot(N.n, S.a) < 0 ? mul(N.n, -1) : N.n;
        const pts = vertsOf(i), m = pts.reduce((s, p) => add(s, mul(p, 1 / pts.length)), [0, 0, 0]);
        // Girar el eje alrededor del punto del eje más cercano a la región.
        const c0 = add(S.c, mul(S.a, dot(sub(m, S.c), S.a)));
        const c1 = add(c0, mul(a, -dot(sub(m, S.c), S.a)));
        const hc = dot(N.n, c1) - N.d;
        if (Math.abs(Math.abs(hc) - S.r) > 3 * tol) continue;
        const sg = Math.sign(hc) || 1;
        const options = [
          { ...S, a, c: add(c1, mul(N.n, sg * S.r - hc)) },
          { ...S, a, c: c1, r: Math.abs(hc) },
        ];
        if (Math.abs(Math.abs(hc) - S.r) < 1e-12 && Math.abs(dot(N.n, S.a)) > 1 - 1e-12) continue;
        for (const S2 of options) if (set(i, S2)) break;
      }
    }
  }

  // 3b) Redondeo entre un plano y un cilindro perpendicular a él (borde de un taladro o de un
  //     saliente): el toro es tangente a ambos. Queda un solo parámetro, el radio del tubo, que
  //     se elige para seguir la malla lo mejor posible.
  for (const i of order) {
    const S = prims[i].S; if (S.type !== 'torus') continue;
    let P = null, Cy = null;
    for (const j of nb[i]) {
      if (!live(j)) continue; const N = prims[j].S;
      if (N.type === 'plane' && Math.abs(dot(N.n, S.a)) > Math.cos(1.5 * RAD) && Math.abs(Math.abs(dot(N.n, S.c) - N.d) - S.r) < 3 * tol) P = N;
      if (N.type === 'cylinder' && Math.abs(dot(N.a, S.a)) > Math.cos(1.5 * RAD) && lineDist(S.c, N.p, N.a) < 3 * tol && Math.min(Math.abs(S.R - N.r - S.r), Math.abs(S.R - N.r + S.r)) < 3 * tol) Cy = N;
    }
    if (!P || !Cy || Math.abs(dot(P.n, Cy.a)) < Math.cos(1.5 * RAD)) continue;
    const a = dot(P.n, S.a) < 0 ? mul(P.n, -1) : P.n;
    if (Math.abs(dot(P.n, Cy.a)) < 1 - 1e-12) {
      // El cilindro también pasa a ser perpendicular al plano.
      const k = prims.findIndex(q => q.S === Cy), pts = vertsOf(k), m = pts.reduce((s2, p) => add(s2, mul(p, 1 / pts.length)), [0, 0, 0]);
      const pivot = add(Cy.p, mul(Cy.a, dot(sub(m, Cy.p), Cy.a)));
      if (!set(k, { ...Cy, a, p: pivot })) continue;
      Cy = prims[k].S;
    }
    const t0 = (P.d - dot(P.n, Cy.p)) / dot(P.n, a), base = add(Cy.p, mul(a, t0));   // eje ∩ plano
    const side = Math.sign(dot(sub(S.c, base), a)) || 1, outer = Math.abs(S.R - Cy.r - S.r) < Math.abs(S.R - Cy.r + S.r);
    const make = r => ({ ...S, a, c: add(base, mul(a, side * r)), r, R: outer ? Cy.r + r : Cy.r - r });
    const pts = vertsOf(i), err = r => { const T2 = make(r); return pts.reduce((m, p) => Math.max(m, Math.abs(dist(T2, p))), 0); };
    let lo = Math.max(1e-3, S.r - 4 * tol), hi = S.r + 4 * tol;
    for (let k = 0; k < 40; k++) { const m1 = lo + (hi - lo) / 3, m2 = hi - (hi - lo) / 3; if (err(m1) < err(m2)) hi = m2; else lo = m1; }
    set(i, make((lo + hi) / 2));
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
  // 4b) Planos vecinos casi coplanarios (< 2°): abanicos de triángulos de la transición entre dos
  //     caras se habían reconocido como planos aparte. Se unen si el plano común encaja.
  for (let pass = 0, again = true; again && pass < 6 && !globalThis.process?.env?.NOCOPLANAR; pass++) { again = false; for (const i of order) {
    if (!live(i) || prims[i].S.type !== 'plane') continue;
    for (const j of nb[i]) {
      if (j === i || !live(j) || prims[j].S.type !== 'plane' || area[j] < area[i]) continue;
      const a = prims[i].S, b = prims[j].S, c = dot(a.n, b.n);
      if (c < Math.cos(2 * RAD)) continue;
      const pts = vertsOf(i).concat(vertsOf(j)), F = fitPlane(pts);
      if (pts.some(p => Math.abs(dot(F.n, p) - F.d) > 2 * tol)) continue;
      const n = dot(F.n, b.n) < 0 ? mul(F.n, -1) : F.n;
      if (globalThis.COPDEBUG) console.log('coplanar', i, j, area[i].toFixed(2), area[j].toFixed(2), (Math.acos(Math.min(1, c)) / RAD).toFixed(3), (Math.acos(Math.min(1, dot(n, b.n))) / RAD).toFixed(3));
      prims[j] = { ...prims[j], S: { type: 'plane', n, d: dot(n, F.o) } };
      for (const t of tris[i]) { reg[t] = j; tris[j].push(t); }
      area[j] += area[i]; area[i] = 0; tris[i] = []; changes++; again = true;
      for (const k of nb[i]) if (k !== j) { nb[j].add(k); nb[k].add(j); }
      break;
    }
  } }

  // 4c) Triángulos degenerados (casi sin anchura) que cayeron en una región cuya superficie es
  //     perpendicular a ellos: están sobre las dos superficies a la vez, pero dejarlos ahí crea
  //     aristas superpuestas que Parasolid rechaza (Onshape llega a eliminar el trozo de sólido).
  if (!globalThis.process?.env?.NOREASSIGN) for (let pass = 0; pass < 3; pass++) {
    let moved = 0;
    for (let t = 0; t < nt; t++) {
      const r = reg[t]; if (!live(r)) continue;
      const P = [0, 1, 2].map(k => { const v = T[t * 3 + k]; return [V[v * 3], V[v * 3 + 1], V[v * 3 + 2]]; });
      const n = cross(sub(P[1], P[0]), sub(P[2], P[0])), l = norm(n); if (!(l > 0)) continue;
      // Solo triángulos degenerados: altura menor que la tolerancia.
      const longest = Math.max(norm(sub(P[1], P[0])), norm(sub(P[2], P[1])), norm(sub(P[0], P[2])));
      if (l / longest > tol) continue;
      const c = mul(add(add(P[0], P[1]), P[2]), 1 / 3), u = mul(n, 1 / l);
      if (Math.abs(dot(u, normalAt(prims[r].S, c))) > 0.5) continue;
      let to = -1;
      for (let k = 0; k < 3; k++) {
        const h = adj[t * 3 + k]; if (h < 0) continue; const q = reg[(h / 3) | 0];
        // Solo entre planos: con superficies curvas el cambio altera contornos ya correctos.
        if (q === r || !live(q) || prims[q].S.type !== 'plane' || prims[r].S.type !== 'plane') continue;
        const S = prims[q].S;
        if (Math.abs(dot(u, normalAt(S, c))) > 0.85 && P.every(p => Math.abs(dist(S, p)) <= tol)) { to = q; break; }
      }
      if (to < 0) continue;
      if (globalThis.REDEBUG) console.log('reasigna', t, prims[r].S.type, '->', prims[to].S.type, (l / longest).toExponential(1), c.map(x => x.toFixed(2)).join(','));
      reg[t] = to; tris[r].splice(tris[r].indexOf(t), 1); tris[to].push(t); area[r] -= mesh.A[t]; area[to] += mesh.A[t]; moved++; changes++;
    }
    if (!moved) break;
  }

  // 5) Franjas planas diminutas junto a un redondeo que en realidad son facetas de él (la fila
  //    de triángulos que toca la cara tangente): se integran en la superficie curva.
  //    Se repite porque al unir uno aparecen vecinos nuevos.
  for (let pass = 0, again = true; again && pass < 6 && !globalThis.process?.env?.NOABSORB; pass++) { again = false; for (const i of order) {
    if (!live(i) || prims[i].S.type !== 'plane') continue;
    let best = null;
    for (const j of nb[i]) {
      if (!live(j) || prims[j].S.type === 'plane' || area[i] > 0.05 * area[j]) continue;
      const w = worst(i, prims[j].S);
      if (!(w <= 1.5 * tol) || (best && w >= best.w)) continue;
      // Cada triángulo debe estar orientado como la superficie curva en ese punto.
      const S = prims[j].S;
      const aligned = tris[i].every(t => {
        const P = [0, 1, 2].map(k => { const v = T[t * 3 + k]; return [V[v * 3], V[v * 3 + 1], V[v * 3 + 2]]; });
        const n = cross(sub(P[1], P[0]), sub(P[2], P[0])), l = norm(n); if (!(l > 0)) return true;
        const c = mul(add(add(P[0], P[1]), P[2]), 1 / 3);
        return Math.abs(dot(mul(n, 1 / l), normalAt(S, c))) > Math.cos(15 * RAD);
      });
      if (aligned) best = { j, w };
    }
    if (!best) continue;
    for (const t of tris[i]) { reg[t] = best.j; tris[best.j].push(t); }
    area[best.j] += area[i]; area[i] = 0; tris[i] = []; changes++; again = true;
    for (const k of nb[i]) if (k !== best.j) { nb[best.j].add(k); nb[k].add(best.j); }
  } }
  return changes;
}
void cross; void unit;
