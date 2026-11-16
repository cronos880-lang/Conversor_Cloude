// Segmentación de la malla en regiones que pertenecen a una misma superficie analítica.
import { dist, normalAt, refine, fitSphere, fitCylinder, fitCone, fitTorus, fitPlane, dot, sub, cross, unit, perp, DOF } from './fit.js';

const RAD = Math.PI / 180;
const TYPES = ['cylinder', 'cone', 'sphere', 'torus'];

export function defaultTolerance(mesh) {
  let maxAbs = 0;
  for (let j = 0; j < 3; j++) maxAbs = Math.max(maxAbs, Math.abs(mesh.min[j]), Math.abs(mesh.max[j]));
  return Math.max(mesh.diag * 1e-5, maxAbs * 2e-6);
}

/**
 * Devuelve { reg, prims }: reg[t] es el índice de región del triángulo t y
 * prims[i] = { S, sign } la superficie de la región i (sign = +1 si la normal
 * exterior de la pieza coincide con la normal natural de la superficie).
 */
export function segment(mesh, options = {}) {
  const { V, T, N, A, adj, nt } = mesh;
  const tol = options.tol ?? defaultTolerance(mesh);
  const cosSmooth = Math.cos((options.smoothAngle ?? 32) * RAD);
  const cosNormal = Math.cos(28 * RAD);
  const reg = new Int32Array(nt).fill(-1), prims = [];
  const P = i => [V[i * 3], V[i * 3 + 1], V[i * 3 + 2]];
  const Nt = t => [N[t * 3], N[t * 3 + 1], N[t * 3 + 2]];
  const centroid = t => {
    const a = T[t * 3] * 3, b = T[t * 3 + 1] * 3, c = T[t * 3 + 2] * 3;
    return [(V[a] + V[b] + V[c]) / 3, (V[a + 1] + V[b + 1] + V[c + 1]) / 3, (V[a + 2] + V[b + 2] + V[c + 2]) / 3];
  };
  const valence = new Int32Array(mesh.nv);
  for (let i = 0; i < nt * 3; i++) valence[T[i]]++;
  const nbr = (t, k) => { const h = adj[t * 3 + k]; return h < 0 ? -1 : (h / 3) | 0; };
  const sliver = t => {
    const a = P(T[t * 3]), b = P(T[t * 3 + 1]), c = P(T[t * 3 + 2]);
    const longest = Math.max(Math.hypot(...sub(a, b)), Math.hypot(...sub(b, c)), Math.hypot(...sub(c, a)));
    return 2 * A[t] < longest * tol * 2;
  };

  // ---- Facetas coplanarias (sobre los triángulos con reg === -1).
  function planeFacets(assign) {
    const facet = new Int32Array(nt).fill(-1), list = [];
    const order = [];
    for (let t = 0; t < nt; t++) if (reg[t] < 0) order.push(t);
    order.sort((a, b) => A[b] - A[a]);
    const stack = [];
    for (const s of order) {
      if (facet[s] >= 0 || !(A[s] > 0)) continue;
      const n = Nt(s), d = dot(n, P(T[s * 3])), id = list.length, tris = [s];
      facet[s] = id; stack.push(s);
      while (stack.length) {
        const t = stack.pop();
        for (let k = 0; k < 3; k++) {
          const u = nbr(t, k); if (u < 0 || reg[u] >= 0 || facet[u] >= 0) continue;
          let ok = true;
          for (let j = 0; j < 3 && ok; j++) { const i = T[u * 3 + j] * 3; ok = Math.abs(n[0] * V[i] + n[1] * V[i + 1] + n[2] * V[i + 2] - d) <= tol; }
          if (!ok) continue;
          if (!sliver(u) && N[u * 3] * n[0] + N[u * 3 + 1] * n[1] + N[u * 3 + 2] * n[2] < 0.94) continue;
          facet[u] = id; tris.push(u); stack.push(u);
        }
      }
      list.push({ n, d, tris });
    }
    // Triángulos sin área que no se hayan absorbido: faceta propia.
    for (const s of order) if (facet[s] < 0) { facet[s] = list.length; list.push({ n: Nt(s), d: dot(Nt(s), P(T[s * 3])), tris: [s] }); }
    if (assign) for (const f of list) {
      const id = prims.length; let S = { type: 'plane', n: f.n, d: f.d };
      if (f.tris.length > 2) {
        const set = new Set(); for (const t of f.tris) for (let j = 0; j < 3; j++) set.add(T[t * 3 + j]);
        const pts = [...set].map(P), fit = fitPlane(pts);
        if (dot(fit.n, f.n) < 0) { fit.n = fit.n.map(x => -x); fit.d = -fit.d; }
        if (dot(fit.n, f.n) > 0.999 && pts.every(p => Math.abs(dot(fit.n, p) - fit.d) <= tol)) S = { type: 'plane', n: fit.n, d: fit.d };
      }
      prims.push({ S, sign: 1 });
      for (const t of f.tris) reg[t] = id;
    }
    return { facet, list };
  }

  // ---- Prueba de pertenencia de un triángulo a una superficie.
  function spanCos(S, t) {
    if (S.type === 'cone') {
      // Cerca del vértice del cono la normal no está definida: medir el giro alrededor del eje.
      const rad = [];
      for (let j = 0; j < 3; j++) {
        const v = sub(P(T[t * 3 + j]), S.apex), h = dot(v, S.a), r = sub(v, S.a.map(x => x * h));
        if (Math.hypot(...r) > tol * 4) rad.push(unit(r));
      }
      let c = 1;
      for (let i = 0; i < rad.length; i++) for (let j = i + 1; j < rad.length; j++) c = Math.min(c, dot(rad[i], rad[j]));
      return c;
    }
    const n0 = normalAt(S, P(T[t * 3])), n1 = normalAt(S, P(T[t * 3 + 1])), n2 = normalAt(S, P(T[t * 3 + 2]));
    return Math.min(dot(n0, n1), dot(n1, n2), dot(n0, n2));
  }
  function fits(S, sign, t, cosSpan) {
    for (let j = 0; j < 3; j++) if (Math.abs(dist(S, P(T[t * 3 + j]))) > tol) return false;
    if (!(A[t] > 0)) return false;
    if (sign * dot(Nt(t), normalAt(S, centroid(t))) < cosNormal) return false;
    return spanCos(S, t) >= cosSpan;
  }
  const mark = new Int32Array(nt); let stamp = 0;
  function grow(S, sign, seeds, cosSpan) {
    stamp++; const out = [], stack = [];
    for (const s of seeds) if (mark[s] !== stamp && reg[s] < 0 && fits(S, sign, s, cosSpan)) { mark[s] = stamp; stack.push(s); out.push(s); break; }
    while (stack.length) {
      const t = stack.pop();
      for (let k = 0; k < 3; k++) {
        const u = nbr(t, k); if (u < 0 || reg[u] >= 0 || mark[u] === stamp) continue;
        mark[u] = stamp;
        if (fits(S, sign, u, cosSpan)) { out.push(u); stack.push(u); }
      }
    }
    return out;
  }
  function vertsOf(tris, limit = Infinity) {
    const set = new Set();
    for (const t of tris) { set.add(T[t * 3]); set.add(T[t * 3 + 1]); set.add(T[t * 3 + 2]); }
    let ids = [...set];
    if (ids.length > limit) { const step = ids.length / limit, pick = []; for (let i = 0; i < limit; i++) pick.push(ids[Math.floor(i * step)]); ids = pick; }
    return ids.map(P);
  }
  const areaOf = tris => { let s = 0; for (const t of tris) s += A[t]; return s; };

  function distinctAngles(values, eps) {
    values.sort((a, b) => a - b); let n = values.length ? 1 : 0;
    for (let i = 1; i < values.length; i++) if (values[i] - values[i - 1] > eps) n++;
    return n;
  }
  /** ¿La región verifica de verdad la superficie (más puntos que grados de libertad)? */
  function verified(S, tris) {
    const pts = vertsOf(tris);
    if (tris.length < 4 || pts.length < DOF[S.type] + 3) return false;
    const big = mesh.diag * 200;
    if (S.type === 'sphere' || S.type === 'torus') {
      // Exigir vértices interiores: una sola banda de triángulos es indistinguible de un cono.
      const count = new Map(); let interior = 0;
      for (const t of tris) for (let j = 0; j < 3; j++) { const v = T[t * 3 + j]; count.set(v, (count.get(v) || 0) + 1); }
      for (const [v, c] of count) if (c === valence[v]) interior++;
      if (interior < 2) return false;
    }
    if (S.type === 'sphere') return S.r < big && S.r > tol * 20;
    const a = S.a, e1 = perp(a), e2 = cross(a, e1), o = S.p || S.apex || S.c;
    const theta = pts.map(p => { const v = sub(p, o); return Math.atan2(dot(v, e2), dot(v, e1)); });
    if (distinctAngles(theta, 1e-3) < 4) return false;
    if (S.type === 'cylinder') return S.r < big && S.r > tol * 20;
    if (S.type === 'cone') return S.alpha > 1 * RAD && S.alpha < 89 * RAD;
    if (S.type === 'torus') {
      if (!(S.r > tol * 20 && S.R > S.r * 1.0001 && S.R < big)) return false;
      const phi = pts.map(p => { const v = sub(p, o), h = dot(v, a); return Math.atan2(h, Math.sqrt(Math.max(0, dot(v, v) - h * h)) - S.R); });
      return distinctAngles(phi, 1e-3) >= 4;
    }
    return false;
  }

  function hypothesis(type, cluster) {
    const pts = vertsOf(cluster), C = cluster.map(centroid), Nrm = cluster.map(Nt);
    let S = type === 'sphere' ? fitSphere(pts) : type === 'cylinder' ? fitCylinder(pts, Nrm) : type === 'cone' ? fitCone(pts, C, Nrm) : fitTorus(pts, C, Nrm);
    if (!S) return null;
    let size = 0; for (const p of pts) size = Math.max(size, Math.hypot(...sub(p, pts[0])));
    S = refine(S, pts, size);
    for (const p of pts) if (!(Math.abs(dist(S, p)) <= tol)) return null;
    let s = 0; for (let i = 0; i < cluster.length; i++) s += A[cluster[i]] * dot(Nrm[i], normalAt(S, C[i]));
    const sign = s >= 0 ? 1 : -1;
    let span = 1; for (const t of cluster) span = Math.min(span, spanCos(S, t));
    const cosSpan = Math.cos(Math.min(35 * RAD, 2.6 * Math.acos(Math.min(1, span)) + 1 * RAD));
    let tris = grow(S, sign, cluster, cosSpan);
    for (let round = 0; round < 4 && tris.length >= 4; round++) {
      const S2 = refine(S, vertsOf(tris, 300), mesh.diag);
      const tris2 = grow(S2, sign, cluster, cosSpan);
      if (tris2.length < tris.length) break;
      const same = tris2.length === tris.length;
      S = S2; tris = tris2;
      if (same) break;
    }
    if (!verified(S, tris)) return null;
    return { S, sign, tris, area: areaOf(tris) };
  }

  function clusterAt(seed, limit, blocked) {
    const out = [seed], seen = new Set([seed]);
    for (let i = 0; i < out.length && out.length < limit; i++) {
      const t = out[i];
      for (let k = 0; k < 3 && out.length < limit; k++) {
        const u = nbr(t, k); if (u < 0 || reg[u] >= 0 || seen.has(u)) continue;
        seen.add(u);
        if (blocked && blocked[u]) continue;
        if (N[t * 3] * N[u * 3] + N[t * 3 + 1] * N[u * 3 + 1] + N[t * 3 + 2] * N[u * 3 + 2] < cosSmooth) continue;
        out.push(u);
      }
    }
    return out;
  }
  const curved = cluster => {
    const n = Nt(cluster[0]);
    return cluster.some(t => dot(n, Nt(t)) < 1 - 1e-7);
  };

  // ---- Pasadas de detección de superficies curvas.
  const { facet, list } = planeFacets(false);
  const blockedBig = new Uint8Array(nt);
  for (let t = 0; t < nt; t++) if (list[facet[t]].tris.length >= 3) blockedBig[t] = 1;
  for (const blocked of [blockedBig, null]) {
    for (let seed = 0; seed < nt; seed++) {
      if (reg[seed] >= 0 || (blocked && blocked[seed]) || !(A[seed] > 0)) continue;
      let best = null;
      const small = clusterAt(seed, 6, blocked), large = clusterAt(seed, 16, blocked);
      if (small.length < 3 || !curved(large)) continue;
      for (const type of TYPES) {
        let h = null;
        if (type === 'cylinder' || type === 'sphere') h = curved(small) ? hypothesis(type, small) : null;
        if (!h && large.length >= 6) h = hypothesis(type, large);
        if (h && (!best || h.area > best.area * 1.25)) best = h;
      }
      if (!best) continue;
      const id = prims.length; prims.push({ S: best.S, sign: best.sign });
      for (const t of best.tris) reg[t] = id;
    }
  }
  planeFacets(true);
  return { reg, prims, tol };
}
