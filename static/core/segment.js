// Segmentación de la malla en regiones que pertenecen a una misma superficie analítica.
import { dist, normalAt, refine, fitSphere, fitCylinder, fitCone, fitConeFromNormals, fitTorus, fitPlane, dot, sub, cross, perp, DOF } from './fit.js';

const RAD = Math.PI / 180;
const TYPES = ['cylinder', 'cone', 'sphere', 'torus'];

/** Tolerancia de pertenencia a una superficie (mm): absorbe el ruido habitual de los STL. */
export function defaultTolerance(mesh) {
  return Math.min(0.03, Math.max(0.01, mesh.diag * 2e-4));
}

/**
 * Devuelve { reg, prims, tol }: reg[t] es el índice de región del triángulo t y
 * prims[i] = { S, sign, res } la superficie de la región i (sign = +1 si la normal
 * exterior de la pieza coincide con la normal natural de la superficie; res es
 * la desviación máxima de los vértices de la región respecto a la superficie).
 */
export function segment(mesh, options = {}) {
  const { V, T, N, A, adj, nt, nv } = mesh;
  const tol = options.tol ?? defaultTolerance(mesh);
  const cosSmooth = Math.cos((options.smoothAngle ?? 32) * RAD);
  const cosNormal = Math.cos(28 * RAD);
  const reg = new Int32Array(nt).fill(-1);
  let prims = [];
  const P = i => [V[i * 3], V[i * 3 + 1], V[i * 3 + 2]];
  const Nt = t => [N[t * 3], N[t * 3 + 1], N[t * 3 + 2]];
  const centroid = t => {
    const a = T[t * 3] * 3, b = T[t * 3 + 1] * 3, c = T[t * 3 + 2] * 3;
    return [(V[a] + V[b] + V[c]) / 3, (V[a + 1] + V[b + 1] + V[c + 1]) / 3, (V[a + 2] + V[b + 2] + V[c + 2]) / 3];
  };
  const valence = new Int32Array(nv);
  for (let i = 0; i < nt * 3; i++) valence[T[i]]++;
  const nbr = (t, k) => { const h = adj[t * 3 + k]; return h < 0 ? -1 : (h / 3) | 0; };
  const longest = new Float64Array(nt);
  for (let t = 0; t < nt; t++) {
    let m = 0;
    for (let k = 0; k < 3; k++) {
      const a = T[t * 3 + k] * 3, b = T[t * 3 + (k + 1) % 3] * 3;
      m = Math.max(m, Math.hypot(V[a] - V[b], V[a + 1] - V[b + 1], V[a + 2] - V[b + 2]));
    }
    longest[t] = m;
  }
  const sliver = (t, eps) => 2 * A[t] < longest[t] * eps * 2;

  // ---- Facetas coplanarias sobre los triángulos sin región.
  function planeFacets(eps, cosN) {
    const facet = new Int32Array(nt).fill(-1), list = [], order = [];
    for (let t = 0; t < nt; t++) if (reg[t] < 0) order.push(t);
    order.sort((a, b) => A[b] - A[a]);
    const stack = [];
    for (const s of order) {
      if (facet[s] >= 0 || !(A[s] > 0)) continue;
      const nx = N[s * 3], ny = N[s * 3 + 1], nz = N[s * 3 + 2], i0 = T[s * 3] * 3;
      const d = nx * V[i0] + ny * V[i0 + 1] + nz * V[i0 + 2], id = list.length, tris = [s];
      facet[s] = id; stack.push(s);
      while (stack.length) {
        const t = stack.pop();
        for (let k = 0; k < 3; k++) {
          const u = nbr(t, k); if (u < 0 || reg[u] >= 0 || facet[u] >= 0) continue;
          let ok = true;
          for (let j = 0; j < 3 && ok; j++) { const i = T[u * 3 + j] * 3; ok = Math.abs(nx * V[i] + ny * V[i + 1] + nz * V[i + 2] - d) <= eps; }
          if (!ok) continue;
          if (!sliver(u, eps) && N[u * 3] * nx + N[u * 3 + 1] * ny + N[u * 3 + 2] * nz < cosN) continue;
          facet[u] = id; tris.push(u); stack.push(u);
        }
      }
      list.push({ n: [nx, ny, nz], d, tris });
    }
    for (const s of order) if (facet[s] < 0) { facet[s] = list.length; list.push({ n: Nt(s), d: dot(Nt(s), P(T[s * 3])), tris: [s] }); }
    return { facet, list };
  }

  // ---- Evaluación (con caché por vértice) de distancia y normal natural a una superficie.
  const vStamp = new Int32Array(nv), vD = new Float64Array(nv), vN = new Float64Array(nv * 3);
  let stamp = 0, cur = null;
  function setSurface(S) { cur = S; stamp++; }
  function evalVertex(i) {
    if (vStamp[i] === stamp) return;
    vStamp[i] = stamp;
    const S = cur, x = V[i * 3], y = V[i * 3 + 1], z = V[i * 3 + 2];
    let d, nx = 0, ny = 0, nz = 0;
    if (S.type === 'sphere') {
      const vx = x - S.c[0], vy = y - S.c[1], vz = z - S.c[2], l = Math.hypot(vx, vy, vz);
      d = l - S.r; if (l > 0) { nx = vx / l; ny = vy / l; nz = vz / l; }
    } else {
      const o = S.p || S.apex || S.c, a = S.a;
      const vx = x - o[0], vy = y - o[1], vz = z - o[2], h = vx * a[0] + vy * a[1] + vz * a[2];
      const rx = vx - h * a[0], ry = vy - h * a[1], rz = vz - h * a[2], rho = Math.hypot(rx, ry, rz);
      const ux = rho > 0 ? rx / rho : 0, uy = rho > 0 ? ry / rho : 0, uz = rho > 0 ? rz / rho : 0;
      if (S.type === 'cylinder') { d = rho - S.r; nx = ux; ny = uy; nz = uz; }
      else if (S.type === 'cone') {
        const ca = Math.cos(S.alpha), sa = Math.sin(S.alpha);
        d = rho * ca - h * sa;
        if (rho > tol * 4) { nx = ux * ca - a[0] * sa; ny = uy * ca - a[1] * sa; nz = uz * ca - a[2] * sa; }
      } else {
        const q = rho - S.R, l = Math.hypot(q, h);
        d = l - S.r;
        if (l > 0) { nx = (ux * q + a[0] * h) / l; ny = (uy * q + a[1] * h) / l; nz = (uz * q + a[2] * h) / l; }
      }
    }
    vD[i] = d; vN[i * 3] = nx; vN[i * 3 + 1] = ny; vN[i * 3 + 2] = nz;
  }
  /** Coseno del giro de la normal de la superficie a lo largo del triángulo. */
  function spanCos(t) {
    const a = T[t * 3], b = T[t * 3 + 1], c = T[t * 3 + 2];
    evalVertex(a); evalVertex(b); evalVertex(c);
    let m = 1;
    const pair = (i, j) => {
      const li = vN[i * 3] * vN[i * 3] + vN[i * 3 + 1] * vN[i * 3 + 1] + vN[i * 3 + 2] * vN[i * 3 + 2];
      const lj = vN[j * 3] * vN[j * 3] + vN[j * 3 + 1] * vN[j * 3 + 1] + vN[j * 3 + 2] * vN[j * 3 + 2];
      if (li < 0.5 || lj < 0.5) return;
      m = Math.min(m, vN[i * 3] * vN[j * 3] + vN[i * 3 + 1] * vN[j * 3 + 1] + vN[i * 3 + 2] * vN[j * 3 + 2]);
    };
    pair(a, b); pair(b, c); pair(a, c);
    return m;
  }
  function fits(sign, t, cosSpan, eps) {
    const a = T[t * 3], b = T[t * 3 + 1], c = T[t * 3 + 2];
    evalVertex(a); if (Math.abs(vD[a]) > eps) return false;
    evalVertex(b); if (Math.abs(vD[b]) > eps) return false;
    evalVertex(c); if (Math.abs(vD[c]) > eps) return false;
    if (!(A[t] > 0)) return false;
    const sx = vN[a * 3] + vN[b * 3] + vN[c * 3], sy = vN[a * 3 + 1] + vN[b * 3 + 1] + vN[c * 3 + 1], sz = vN[a * 3 + 2] + vN[b * 3 + 2] + vN[c * 3 + 2];
    const l = Math.hypot(sx, sy, sz);
    if (!(l > 0)) return false;
    if (!sliver(t, eps) && sign * (N[t * 3] * sx + N[t * 3 + 1] * sy + N[t * 3 + 2] * sz) < cosNormal * l) return false;
    return spanCos(t) >= cosSpan;
  }

  // Facetas planas de protección: una faceta grande solo entra entera en una superficie curva.
  let guard = planeFacets(tol / 3, Math.cos(6 * RAD));
  const mark = new Int32Array(nt); let fMark = new Int32Array(guard.list.length), gStamp = 0;
  /**
   * Crecimiento de una región desde `seeds`. `allow(u)` decide si el triángulo u
   * (posiblemente ya asignado) puede incorporarse.
   */
  function grow(S, sign, seeds, cosSpan, eps, allow) {
    setSurface(S); gStamp++;
    const out = [], stack = [];
    const take = u => {
      const f = guard.facet[u], ft = f >= 0 ? guard.list[f].tris : null;
      if (ft && ft.length > 2) {
        if (fMark[f] === gStamp) return;
        fMark[f] = gStamp;
        const members = ft.filter(x => allow(x));
        for (const x of members) mark[x] = gStamp;
        if (members.length !== ft.length || !members.every(x => fits(sign, x, cosSpan, eps))) return;
        for (const x of members) { out.push(x); stack.push(x); }
      } else {
        mark[u] = gStamp;
        if (fits(sign, u, cosSpan, eps)) { out.push(u); stack.push(u); }
      }
    };
    for (const s of seeds) { if (mark[s] !== gStamp && allow(s)) take(s); if (out.length) break; }
    while (stack.length) {
      const t = stack.pop();
      for (let k = 0; k < 3; k++) {
        const u = nbr(t, k); if (u < 0 || mark[u] === gStamp || !allow(u)) continue;
        take(u);
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
  function distinct(values, eps) {
    values.sort((a, b) => a - b); let n = values.length ? 1 : 0, last = values[0];
    for (let i = 1; i < values.length; i++) if (values[i] - last > eps) { n++; last = values[i]; }
    return n;
  }
  /** ¿La región verifica de verdad la superficie (más datos que grados de libertad)? */
  function verified(S, tris) {
    if (tris.length < 4) return false;
    const pts = vertsOf(tris);
    if (pts.length < DOF[S.type] + 3) return false;
    const big = mesh.diag * 100;
    // Una región plana dentro de la tolerancia no es una superficie curva.
    const flat = fitPlane(pts);
    if (pts.every(p => Math.abs(dot(flat.n, p) - flat.d) <= tol * 1.5)) return false;
    if (S.type === 'sphere' || S.type === 'torus') {
      // Exigir vértices interiores: una sola banda de triángulos es indistinguible de un cono.
      const count = new Map(); let interior = 0;
      for (const t of tris) for (let j = 0; j < 3; j++) { const v = T[t * 3 + j]; count.set(v, (count.get(v) || 0) + 1); }
      for (const [v, c] of count) if (c === valence[v]) interior++;
      if (interior < 2) return false;
    }
    // La curvatura debe apreciarse: flecha de la región varias veces mayor que la tolerancia.
    let extent = 0; for (const p of pts) extent = Math.max(extent, Math.hypot(p[0] - pts[0][0], p[1] - pts[0][1], p[2] - pts[0][2]));
    const sag = r => { const half = Math.min(extent / 2, r); return r - Math.sqrt(r * r - half * half); };
    if (S.type === 'sphere') return S.r < big && S.r > tol * 5 && sag(S.r) > 4 * tol;
    if (S.type === 'cylinder' && !(sag(S.r) > 4 * tol)) return false;
    const a = S.a, e1 = perp(a), e2 = cross(a, e1), o = S.p || S.apex || S.c;
    const rmin = S.type === 'cylinder' ? S.r : S.type === 'torus' ? S.r : 0;
    const theta = pts.map(p => { const v = sub(p, o); return Math.atan2(dot(v, e2), dot(v, e1)); });
    // Separación angular mínima equivalente a varias veces la tolerancia sobre el radio.
    const epsA = Math.max(1e-3, rmin > 0 ? 4 * tol / (S.type === 'torus' ? S.R : rmin) : 1e-3);
    if (distinct(theta, epsA) < 4) return false;
    if (S.type === 'cylinder') return S.r < big && S.r > tol * 5;
    if (S.type === 'cone') return S.alpha > 1 * RAD && S.alpha < 89 * RAD;
    if (!(S.r > tol * 5 && S.R > S.r * 1.0001 && S.R < big)) return false;
    const phi = pts.map(p => { const v = sub(p, o), h = dot(v, a); return Math.atan2(h, Math.sqrt(Math.max(0, dot(v, v) - h * h)) - S.R); });
    return distinct(phi, Math.max(1e-3, 4 * tol / S.r)) >= 4;
  }
  const residual = (S, pts) => { let m = 0; for (const p of pts) m = Math.max(m, Math.abs(dist(S, p))); return m; };
  const free = u => reg[u] < 0;

  function hypothesis(type, cluster) {
    const pts = vertsOf(cluster), C = cluster.map(centroid), Nrm = cluster.map(Nt);
    let size = 0; for (const p of pts) size = Math.max(size, Math.hypot(...sub(p, pts[0])));
    const inits = type === 'sphere' ? [fitSphere(pts)] : type === 'cylinder' ? [fitCylinder(pts, Nrm)]
      : type === 'cone' ? [fitCone(pts, C, Nrm), fitConeFromNormals(pts, Nrm, C)] : [fitTorus(pts, C, Nrm)];
    let S = null, best = Infinity;
    for (const init of inits) {
      if (!init) continue;
      const cand = refine(init, pts, size), r = residual(cand, pts);
      if (r < best) { best = r; S = cand; }
    }
    if (!S) return null;
    // Si el grupo mezcla dos superficies, descartar los triángulos peor ajustados y reintentar.
    for (let round = 0; round < 3 && best > tol; round++) {
      const worst = cluster.map(t => Math.max(...[0, 1, 2].map(j => Math.abs(dist(S, P(T[t * 3 + j]))))));
      const cut = worst.slice().sort((a, b) => a - b)[Math.floor(worst.length * 0.6)];
      const keep = cluster.filter((t, i) => worst[i] <= Math.max(cut, tol));
      if (keep.length === cluster.length || keep.length < 4) return null;
      return hypothesisOn(type, keep, S);
    }
    return finish(type, cluster, S, pts, C, Nrm);
  }
  function hypothesisOn(type, cluster, S0) {
    const pts = vertsOf(cluster), C = cluster.map(centroid), Nrm = cluster.map(Nt);
    if (pts.length < DOF[type] + 1) return null;
    let size = 0; for (const p of pts) size = Math.max(size, Math.hypot(...sub(p, pts[0])));
    let S = refine(S0, pts, size), best = residual(S, pts);
    for (let round = 0; round < 2 && best > tol; round++) {
      const worst = cluster.map(t => Math.max(...[0, 1, 2].map(j => Math.abs(dist(S, P(T[t * 3 + j]))))));
      const cut = worst.slice().sort((a, b) => a - b)[Math.floor(worst.length * 0.6)];
      const keep = cluster.filter((t, i) => worst[i] <= Math.max(cut, tol));
      if (keep.length === cluster.length || keep.length < 4) return null;
      cluster = keep;
      const p2 = vertsOf(cluster); if (p2.length < DOF[type] + 1) return null;
      S = refine(S, p2, size); best = residual(S, p2);
    }
    if (best > tol) return null;
    return finish(type, cluster, S, vertsOf(cluster), cluster.map(centroid), cluster.map(Nt));
  }
  function finish(type, cluster, S, pts, C, Nrm) {
    let s = 0; for (let i = 0; i < cluster.length; i++) s += A[cluster[i]] * dot(Nrm[i], normalAt(S, C[i]));
    const sign = s >= 0 ? 1 : -1;
    setSurface(S);
    let span = 1; for (const t of cluster) span = Math.min(span, spanCos(t));
    const cosSpan = Math.cos(Math.min(35 * RAD, 2.6 * Math.acos(Math.min(1, span)) + 1 * RAD));
    let tris = grow(S, sign, cluster, cosSpan, tol, free);
    for (let round = 0; round < 4 && tris.length >= 4; round++) {
      const S2 = refine(S, vertsOf(tris, 300), mesh.diag);
      const tris2 = grow(S2, sign, cluster, cosSpan, tol, free);
      if (tris2.length < tris.length) break;
      const same = tris2.length === tris.length;
      S = S2; tris = tris2;
      if (same) break;
    }
    if (options.debug) console.log('   ', type, 'cluster', cluster.length, 'grown', tris.length, 'verified', verified(S, tris), S.r && S.r.toFixed(4), 'span', (Math.acos(cosSpan) / RAD).toFixed(1));
    if (!verified(S, tris)) return null;
    return { S, sign, tris, area: areaOf(tris), cosSpan };
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
  // El grupo es curvo si algún vértice se separa del plano del primer triángulo.
  const curved = cluster => {
    const s = cluster[0], n = Nt(s), d = dot(n, P(T[s * 3]));
    for (const t of cluster) for (let j = 0; j < 3; j++) if (Math.abs(dot(n, P(T[t * 3 + j])) - d) > tol) return true;
    return false;
  };

  // ---- Detección de superficies curvas.
  const blockedBig = new Uint8Array(nt);
  for (let t = 0; t < nt; t++) if (guard.list[guard.facet[t]].tris.length >= 7) blockedBig[t] = 1;
  const regions = []; // { S, sign, tris, area, cosSpan }
  const tried = new Uint8Array(nt);
  function seedPass(blocked) {
    for (let seed = 0; seed < nt; seed++) {
      if (reg[seed] >= 0 || tried[seed] === (blocked ? 1 : 2) || (blocked && blocked[seed]) || !(A[seed] > 0)) continue;
      tried[seed] = blocked ? 1 : 2;
      const small = clusterAt(seed, 6, blocked), large = clusterAt(seed, 16, blocked);
      if (options.debug && options.debug(centroid(seed))) console.log('seed', seed, blocked ? 'pass1' : 'pass2', 'small', small.length, 'large', large.length, 'curved', curved(large)); else if (options.debug) continue;
      if (small.length < 3 || !curved(large)) continue;
      let best = null;
      for (const type of TYPES) {
        let h = null;
        if (type === 'cylinder' || type === 'sphere') h = curved(small) ? hypothesis(type, small) : null;
        if (!h && large.length >= 6) h = hypothesis(type, large);
        if (h && (!best || h.area > best.area * 1.25)) best = h;
      }
      if (!best) continue;
      const id = regions.length; regions.push(best);
      for (const t of best.tris) reg[t] = id;
    }
  }
  seedPass(blockedBig);
  seedPass(null);

  // ---- Fusión: las regiones grandes absorben triángulos de regiones curvas menores que encajan.
  const order = regions.map((_, i) => i).sort((a, b) => regions[b].area - regions[a].area);
  let changed = false;
  for (const id of order) {
    const R = regions[id]; if (!R.tris.length) continue;
    const allow = u => reg[u] < 0 || reg[u] === id || (regions[reg[u]].area < R.area);
    for (const t of R.tris) guard.facet[t] = -1;
    const tris = grow(R.S, R.sign, R.tris, R.cosSpan, tol, allow);
    if (tris.length <= R.tris.length) continue;
    const victims = new Set();
    for (const t of tris) { if (reg[t] >= 0 && reg[t] !== id) victims.add(reg[t]); reg[t] = id; }
    const kept = new Set(tris);
    for (const t of R.tris) if (!kept.has(t)) reg[t] = -1;
    R.tris = tris; R.area = areaOf(tris);
    const S2 = refine(R.S, vertsOf(tris, 400), mesh.diag);
    if (residual(S2, vertsOf(tris)) <= tol) R.S = S2;
    for (const v of victims) {
      const W = regions[v], rest = W.tris.filter(t => reg[t] === v);
      // Lo que quede de la región absorbida se libera para volver a analizarlo.
      for (const t of rest) reg[t] = -1;
      W.tris = []; W.area = 0; changed = true;
    }
  }
  if (changed) { tried.fill(0); guard = planeFacets(tol / 3, Math.cos(6 * RAD)); fMark = new Int32Array(guard.list.length); seedPass(null); }

  // ---- Numeración final: regiones curvas y después las facetas planas del resto.
  reg.fill(-1);
  for (const R of regions) {
    if (!R.tris.length) continue;
    const id = prims.length;
    prims.push({ S: R.S, sign: R.sign, res: residual(R.S, vertsOf(R.tris)) });
    for (const t of R.tris) reg[t] = id;
  }
  const { list } = planeFacets(tol, 0.94);
  for (const f of list) {
    const id = prims.length; let S = { type: 'plane', n: f.n, d: f.d };
    const set = new Set(); for (const t of f.tris) for (let j = 0; j < 3; j++) set.add(T[t * 3 + j]);
    const pts = [...set].map(P);
    if (f.tris.length > 1) {
      const fit = fitPlane(pts);
      if (dot(fit.n, f.n) < 0) { fit.n = fit.n.map(x => -x); fit.d = -fit.d; }
      if (dot(fit.n, f.n) > 0.99 && residual(fit, pts) <= residual(S, pts)) S = { type: 'plane', n: fit.n, d: fit.d };
    }
    prims.push({ S, sign: 1, res: residual(S, pts) });
    for (const t of f.tris) reg[t] = id;
  }
  absorbSlivers();
  return { reg, prims, tol };

  /** Las regiones formadas solo por triángulos sin anchura se unen a la región vecina. */
  function absorbSlivers() {
    const tris = prims.map(() => []);
    for (let t = 0; t < nt; t++) tris[reg[t]].push(t);
    for (let r = 0; r < prims.length; r++) {
      if (!tris[r].length || !tris[r].every(t => sliver(t, tol))) continue;
      const share = new Map();
      for (const t of tris[r]) for (let k = 0; k < 3; k++) {
        const u = nbr(t, k); if (u < 0 || reg[u] === r) continue;
        const a = T[t * 3 + k] * 3, b = T[t * 3 + (k + 1) % 3] * 3;
        share.set(reg[u], (share.get(reg[u]) || 0) + Math.hypot(V[a] - V[b], V[a + 1] - V[b + 1], V[a + 2] - V[b + 2]));
      }
      let best = -1, len = 0;
      for (const [q, l] of share) if (l > len) { len = l; best = q; }
      if (best < 0) continue;
      for (const t of tris[r]) { reg[t] = best; tris[best].push(t); }
      tris[r] = [];
    }
  }
}

/**
 * Sustituye las regiones indicadas por facetas planas (una por grupo de
 * triángulos coplanarios). Modifica `seg` y devuelve cuántas regiones cambió.
 */
export function demote(mesh, seg, regions) {
  const { V, T, N, A, adj, nt } = mesh, { reg, prims } = seg;
  const eps = Math.max(mesh.diag * 1e-6, 1e-7);
  const pending = new Uint8Array(nt); let n = 0;
  for (let t = 0; t < nt; t++) if (regions.has(reg[t])) { pending[t] = 1; n++; }
  if (!n) return 0;
  const order = []; for (let t = 0; t < nt; t++) if (pending[t]) order.push(t);
  order.sort((a, b) => A[b] - A[a]);
  const stack = [];
  for (const s of order) {
    if (!pending[s]) continue;
    const nx = N[s * 3], ny = N[s * 3 + 1], nz = N[s * 3 + 2], i0 = T[s * 3] * 3;
    const d = nx * V[i0] + ny * V[i0 + 1] + nz * V[i0 + 2], id = prims.length;
    prims.push({ S: { type: 'plane', n: [nx, ny, nz], d }, sign: 1, res: 0, facet: true });
    pending[s] = 0; reg[s] = id; stack.push(s);
    while (stack.length) {
      const t = stack.pop();
      for (let k = 0; k < 3; k++) {
        const h = adj[t * 3 + k]; if (h < 0) continue;
        const u = (h / 3) | 0; if (!pending[u]) continue;
        let ok = N[u * 3] * nx + N[u * 3 + 1] * ny + N[u * 3 + 2] * nz > 0.9999;
        for (let j = 0; j < 3 && ok; j++) { const i = T[u * 3 + j] * 3; ok = Math.abs(nx * V[i] + ny * V[i + 1] + nz * V[i + 2] - d) <= eps; }
        if (!ok) continue;
        pending[u] = 0; reg[u] = id; stack.push(u);
      }
    }
  }
  return regions.size;
}
