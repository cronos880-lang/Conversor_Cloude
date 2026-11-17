// Superficies extruidas: cadenas de franjas planas, paralelas a una dirección, que juntas
// forman un perfil curvo suave (flancos de engranajes, perfiles de levas, contornos libres).
import { dot, cross, unit, perp, normalAt, dist, fitCircle2D } from './fit.js';

const RAD = Math.PI / 180;

/** Tangentes de una poligonal 2D (dirección de la circunferencia por tres puntos). */
function tangents2(P) {
  const n = P.length, out = [];
  for (let i = 0; i < n; i++) {
    const prev = P[Math.max(0, i - 1)], next = P[Math.min(n - 1, i + 1)];
    if (i > 0 && i < n - 1) {
      const ax = P[i][0] - prev[0], ay = P[i][1] - prev[1], bx = next[0] - P[i][0], by = next[1] - P[i][1];
      const la = Math.hypot(ax, ay), lb = Math.hypot(bx, by);
      const tx = ax * lb / la + bx * la / lb, ty = ay * lb / la + by * la / lb, l = Math.hypot(tx, ty);
      out.push([tx / l, ty / l]);
    } else out.push(null);
  }
  // Extremos: reflejar la tangente vecina respecto a la cuerda.
  for (const [i, j, a, b] of [[0, 1, 0, 1], [n - 1, n - 2, n - 2, n - 1]]) {
    let cx = P[b][0] - P[a][0], cy = P[b][1] - P[a][1]; const l = Math.hypot(cx, cy); cx /= l; cy /= l;
    const t = out[j] || [cx, cy], k = 2 * (t[0] * cx + t[1] * cy);
    out[i] = [k * cx - t[0], k * cy - t[1]];
  }
  return out;
}

/** Spline cúbica C1 (tramos de Bézier) que interpola los puntos 2D. */
function spline2(P) {
  const T = tangents2(P), ctrl = [P[0]], knots = [0], mults = [4];
  let u = 0;
  for (let i = 0; i + 1 < P.length; i++) {
    const cx = P[i + 1][0] - P[i][0], cy = P[i + 1][1] - P[i][1], c = Math.hypot(cx, cy);
    const c0 = (T[i][0] * cx + T[i][1] * cy) / c, c1 = (T[i + 1][0] * cx + T[i + 1][1] * cy) / c;
    const h0 = (2 / 3) * c / (1 + Math.max(0.2, c0)), h1 = (2 / 3) * c / (1 + Math.max(0.2, c1));
    ctrl.push([P[i][0] + T[i][0] * h0, P[i][1] + T[i][1] * h0], [P[i + 1][0] - T[i + 1][0] * h1, P[i + 1][1] - T[i + 1][1] * h1], P[i + 1]);
    u += c; knots.push(u); mults.push(3);
  }
  mults[mults.length - 1] = 4;
  // Muestreo denso para medir distancias.
  const samples = [];
  for (let i = 0; i + 1 < P.length; i++) {
    const [a, b, c, d] = ctrl.slice(i * 3, i * 3 + 4);
    for (let k = i ? 1 : 0; k <= 24; k++) {
      const t = k / 24, s = 1 - t, w0 = s * s * s, w1 = 3 * s * s * t, w2 = 3 * s * t * t, w3 = t * t * t;
      samples.push(w0 * a[0] + w1 * b[0] + w2 * c[0] + w3 * d[0], w0 * a[1] + w1 * b[1] + w2 * c[1] + w3 * d[1]);
    }
  }
  return { ctrl, knots, mults, samples: Float64Array.from(samples) };
}

/**
 * Busca superficies extruidas y sustituye las franjas planas que las forman por una
 * sola región. Modifica `seg` y devuelve cuántas superficies creó.
 */
export function findExtrusions(mesh, seg) {
  const { V, T, N, A, adj, nt } = mesh, { reg, prims, tol } = seg;
  const R = prims.length, area = new Float64Array(R), tris = prims.map(() => []);
  for (let t = 0; t < nt; t++) { area[reg[t]] += A[t]; tris[reg[t]].push(t); }
  // Direcciones candidatas: normales de los planos grandes y ejes de revolución.
  const dirs = [];
  const addDir = d => { if (dirs.length < 12 && !dirs.some(e => Math.abs(dot(e, d)) > 0.9999)) dirs.push(unit(d)); };
  const byArea = prims.map((p, i) => i).filter(i => tris[i].length).sort((a, b) => area[b] - area[a]);
  for (const i of byArea) { const S = prims[i].S; if (S.type === 'plane' && area[i] > 0.01 * mesh.area) addDir(S.n); }
  for (const i of byArea) { const S = prims[i].S; if (S.type === 'cylinder' || S.type === 'cone') addDir(S.a); }

  const used = new Uint8Array(R); let made = 0;
  const same = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1]) <= tol;
  for (const d of dirs) {
    const e1 = perp(d), e2 = cross(d, e1);
    const to2 = i => { const x = V[i * 3], y = V[i * 3 + 1], z = V[i * 3 + 2]; return [x * e1[0] + y * e1[1] + z * e1[2], x * e2[0] + y * e2[1] + z * e2[2]]; };
    // Franjas: planos paralelos a d; sus vértices proyectados forman una poligonal recta ordenada.
    const strip = new Map();
    for (const r of byArea) {
      const S = prims[r].S;
      if (used[r] || S.type !== 'plane' || Math.abs(dot(S.n, d)) > 2e-3 || area[r] > 0.05 * mesh.area) continue;
      const tx = -dot(S.n, e2), ty = dot(S.n, e1); // dirección de la franja en 2D (perpendicular a la normal)
      const pts = [];
      for (const t of tris[r]) for (let j = 0; j < 3; j++) { const p = to2(T[t * 3 + j]); if (!pts.some(q => same(p, q))) pts.push(p); }
      if (pts.length < 2 || pts.length > 64) continue;
      pts.sort((p, q) => (p[0] * tx + p[1] * ty) - (q[0] * tx + q[1] * ty));
      strip.set(r, pts);
    }
    if (strip.size < 3) continue;
    // Enlaces suaves entre franjas vecinas que comparten un extremo.
    const link = new Map(); for (const r of strip.keys()) link.set(r, new Set());
    for (let h = 0; h < nt * 3; h++) {
      const a = reg[(h / 3) | 0], b = adj[h] < 0 ? -1 : reg[(adj[h] / 3) | 0];
      if (a >= b || !strip.has(a) || !strip.has(b)) continue;
      if (dot(prims[a].S.n, prims[b].S.n) < Math.cos(35 * RAD)) continue;
      const pa = strip.get(a), pb = strip.get(b);
      const ends = [pa[0], pa[pa.length - 1]], shared = ends.some(p => same(p, pb[0]) || same(p, pb[pb.length - 1]));
      if (shared) { link.get(a).add(b); link.get(b).add(a); }
    }
    for (const [r, s] of link) if (s.size > 2) { for (const q of s) link.get(q).delete(r); s.clear(); }
    const seen = new Set();
    for (const start of strip.keys()) {
      if (seen.has(start) || link.get(start).size !== 1) continue;
      // Recorrer la cadena desde un extremo.
      const chain = [start]; seen.add(start);
      for (let cur = start; ;) {
        const next = [...link.get(cur)].find(q => !seen.has(q)); if (next === undefined) break;
        chain.push(next); seen.add(next); cur = next;
      }
      if (chain.length < 3) continue;
      // Perfil: puntos de cada franja en el orden de la cadena.
      let profile = null, ok = true;
      for (let k = 0; k < chain.length && ok; k++) {
        let pts = strip.get(chain[k]);
        if (!profile) {
          const nb = strip.get(chain[1]);
          if (!(same(pts[pts.length - 1], nb[0]) || same(pts[pts.length - 1], nb[nb.length - 1]))) pts = pts.slice().reverse();
          profile = pts.slice();
        } else {
          const last = profile[profile.length - 1];
          if (same(pts[pts.length - 1], last)) pts = pts.slice().reverse();
          if (!same(pts[0], last)) { ok = false; break; }
          profile.push(...pts.slice(1));
        }
      }
      if (!ok || profile.length < 4) continue;
      // Suavidad del perfil y giro total apreciable (si no, es un plano).
      let total = 0;
      for (let i = 1; i + 1 < profile.length; i++) {
        const ax = profile[i][0] - profile[i - 1][0], ay = profile[i][1] - profile[i - 1][1], bx = profile[i + 1][0] - profile[i][0], by = profile[i + 1][1] - profile[i][1];
        const ang = Math.abs(Math.atan2(ax * by - ay * bx, ax * bx + ay * by));
        if (ang > 35 * RAD) { ok = false; break; }
        total += ang;
      }
      if (!ok || total < 10 * RAD) continue;
      // Perfil en arco de circunferencia: es un cilindro, no una extrusión libre.
      const circ = fitCircle2D(profile);
      if (circ && profile.every(p => Math.abs(Math.hypot(p[0] - circ[0], p[1] - circ[1]) - circ[2]) <= tol)) {
        const c3 = [e1[0] * circ[0] + e2[0] * circ[1], e1[1] * circ[0] + e2[1] * circ[1], e1[2] * circ[0] + e2[2] * circ[1]];
        const C = { type: 'cylinder', p: c3, a: d, r: circ[2] };
        let s2 = 0;
        for (const r of chain) for (const t of tris[r]) { const c = [0, 1, 2].map(k => (V[T[t * 3] * 3 + k] + V[T[t * 3 + 1] * 3 + k] + V[T[t * 3 + 2] * 3 + k]) / 3); s2 += A[t] * dot([N[t * 3], N[t * 3 + 1], N[t * 3 + 2]], normalAt(C, c)); }
        const r0 = chain[0]; prims[r0] = { S: C, sign: s2 >= 0 ? 1 : -1, res: 0 };
        for (const r of chain) { used[r] = 1; if (r !== r0) for (const t of tris[r]) reg[t] = r0; }
        made++; continue;
      }
      const sp = spline2(profile);
      // La curva no puede separarse de las franjas (cuerdas del perfil): si no, el perfil
      // no es una curva suave sino un polígono con caras planas de verdad.
      { let far = 0;
        for (let i = 0; i + 1 < profile.length; i++) {
          const mx = (profile[i][0] + profile[i + 1][0]) / 2, my = (profile[i][1] + profile[i + 1][1]) / 2, q = sp.samples, k = (i * 24 + 12) * 2;
          far = Math.max(far, Math.hypot(q[k] - mx, q[k + 1] - my));
        }
        if (far > 6 * tol) continue; }
      let h0 = Infinity;
      for (const r of chain) for (const t of tris[r]) for (let j = 0; j < 3; j++) { const i = T[t * 3 + j]; h0 = Math.min(h0, V[i * 3] * d[0] + V[i * 3 + 1] * d[1] + V[i * 3 + 2] * d[2]); }
      const o = [d[0] * h0, d[1] * h0, d[2] * h0], P3 = c => [o[0] + e1[0] * c[0] + e2[0] * c[1], o[1] + e1[1] * c[0] + e2[1] * c[1], o[2] + e1[2] * c[0] + e2[2] * c[1]];
      const S = { type: 'extrusion', d, e1, e2, o, samples: sp.samples, curve: { type: 'bspline', degree: 3, ctrl: sp.ctrl.map(P3), knots: sp.knots, mults: sp.mults } };
      // Todos los vértices deben quedar sobre la superficie.
      let res = 0, s = 0;
      for (const r of chain) for (const t of tris[r]) {
        for (let j = 0; j < 3; j++) { const i = T[t * 3 + j]; res = Math.max(res, Math.abs(dist(S, [V[i * 3], V[i * 3 + 1], V[i * 3 + 2]]))); }
        const c = [0, 1, 2].map(k => (V[T[t * 3] * 3 + k] + V[T[t * 3 + 1] * 3 + k] + V[T[t * 3 + 2] * 3 + k]) / 3);
        s += A[t] * dot([N[t * 3], N[t * 3 + 1], N[t * 3 + 2]], normalAt(S, c));
      }
      if (!(res <= tol)) continue;
      const r0 = chain[0];
      prims[r0] = { S, sign: s >= 0 ? 1 : -1, res };
      for (const r of chain) { used[r] = 1; if (r !== r0) for (const t of tris[r]) reg[t] = r0; }
      made++;
    }
  }
  return made;
}
