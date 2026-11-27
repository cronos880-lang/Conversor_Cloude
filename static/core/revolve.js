// Superficies de revolución coaxiales vistas como perfiles en el semiplano (r, z) de su eje común:
// cilindros, conos y planos perpendiculares son rectas; toros y esferas, circunferencias. Sus
// intersecciones son circunferencias exactas en 3D y las tangencias de los redondeos se imponen en 2D.
import { dot, sub, add, mul, norm } from './fit.js';

const isRev = S => S.type === 'cylinder' || S.type === 'cone' || S.type === 'torus' || S.type === 'sphere';
const axisPoint = S => S.p || S.apex || S.c;

/** Eje común (o, a) de dos superficies, o null si no son coaxiales. */
export function commonAxis(A, B, eps = 1e-6) {
  const R = isRev(A) && A.type !== 'sphere' ? A : isRev(B) && B.type !== 'sphere' ? B : null;
  if (!R) return null;
  const o = axisPoint(R), a = R.a;
  for (const S of [A, B]) {
    if (S === R) continue;
    if (S.type === 'plane') { if (Math.abs(Math.abs(dot(S.n, a)) - 1) > 1e-9) return null; continue; }
    if (S.type === 'sphere') { const v = sub(S.c, o); if (norm(sub(v, mul(a, dot(v, a)))) > eps * (1 + 0.1 * norm(v))) return null; continue; }
    if (!isRev(S) || Math.abs(Math.abs(dot(S.a, a)) - 1) > 1e-9) return null;
    const v = sub(axisPoint(S), o); if (norm(sub(v, mul(a, dot(v, a)))) > eps * (1 + 0.1 * norm(v))) return null;
  }
  return { o, a };
}

/** Perfil de S en el plano (r, z) del eje (o, a). */
export function profileOf(S, o, a) {
  const z = p => dot(sub(p, o), a);
  switch (S.type) {
    case 'cylinder': return { kind: 'line', q: [S.r, 0], t: [0, 1] };
    case 'plane': return { kind: 'line', q: [0, (S.d - dot(S.n, o)) / dot(S.n, a)], t: [1, 0] };
    case 'cone': { const s = Math.sign(dot(S.a, a)) || 1; return { kind: 'line', q: [0, z(S.apex)], t: [Math.sin(S.alpha), s * Math.cos(S.alpha)] }; }
    case 'torus': return { kind: 'circle', c: [S.R, z(S.c)], r: S.r };
    case 'sphere': return { kind: 'circle', c: [0, z(S.c)], r: S.r };
  }
  return null;
}

const d2 = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1]);
const lineNormal = L => [-L.t[1], L.t[0]];

/** Puntos (r, z) comunes a dos perfiles, con `tangent` (radio de curvatura del contacto) si se tocan. */
export function profileMeet(P, Q) {
  if (P.kind === 'circle' && Q.kind === 'line') [P, Q] = [Q, P];
  const out = [];
  if (P.kind === 'line' && Q.kind === 'line') {
    const det = P.t[0] * Q.t[1] - P.t[1] * Q.t[0]; if (Math.abs(det) < 1e-12) return out;
    const w = [Q.q[0] - P.q[0], Q.q[1] - P.q[1]], s = (w[0] * Q.t[1] - w[1] * Q.t[0]) / det;
    out.push({ p: [P.q[0] + s * P.t[0], P.q[1] + s * P.t[1]] });
  } else if (P.kind === 'line' && Q.kind === 'circle') {
    const n = lineNormal(P), h = (Q.c[0] - P.q[0]) * n[0] + (Q.c[1] - P.q[1]) * n[1];
    const foot = [Q.c[0] - n[0] * h, Q.c[1] - n[1] * h], w2 = Q.r * Q.r - h * h;
    if (Math.abs(Math.abs(h) - Q.r) < 1e-9 * Math.max(1, Q.r)) out.push({ p: foot, tangent: Q.r });
    else if (w2 > 0) { const w = Math.sqrt(w2); out.push({ p: [foot[0] + P.t[0] * w, foot[1] + P.t[1] * w] }, { p: [foot[0] - P.t[0] * w, foot[1] - P.t[1] * w] }); }
  } else if (P.kind === 'circle' && Q.kind === 'circle') {
    const d = d2(P.c, Q.c); if (!(d > 0)) return out;
    const x = (d * d + P.r * P.r - Q.r * Q.r) / (2 * d), y2 = P.r * P.r - x * x;
    const u = [(Q.c[0] - P.c[0]) / d, (Q.c[1] - P.c[1]) / d], v = [-u[1], u[0]], m = [P.c[0] + u[0] * x, P.c[1] + u[1] * x];
    if (Math.abs(y2) < 1e-9 * Math.max(1, P.r * P.r)) out.push({ p: m, tangent: Math.min(P.r, Q.r) });
    else if (y2 > 0) { const y = Math.sqrt(y2); out.push({ p: [m[0] + v[0] * y, m[1] + v[1] * y] }, { p: [m[0] - v[0] * y, m[1] - v[1] * y] }); }
  }
  return out.filter(k => k.p[0] > 1e-9);
}

/** Circunferencias 3D { c, a, r, tangent? } donde se cortan dos superficies coaxiales. */
export function coaxialCircles(A, B) {
  const ax = commonAxis(A, B); if (!ax) return [];
  const PA = profileOf(A, ax.o, ax.a), PB = profileOf(B, ax.o, ax.a); if (!PA || !PB) return [];
  return profileMeet(PA, PB).map(k => ({ c: add(ax.o, mul(ax.a, k.p[1])), a: ax.a, r: k.p[0], tangent: k.tangent }));
}

/**
 * Redondeo (toro) tangente a uno o dos perfiles rectos coaxiales (cono, cilindro o plano):
 * devuelve la función radio -> toro, o null. Con dos vecinos el centro queda en el cruce de sus
 * paralelas a distancia ρ; con uno, el centro se desliza por la normal hasta distar exactamente r.
 */
export function filletMaker(T, neighbours) {
  const lines = [];
  for (const N of neighbours) {
    const ax = commonAxis(T, N); if (!ax) continue;
    const L = profileOf(N, T.c, T.a); if (!L || L.kind !== 'line') continue;
    lines.push(L);
  }
  if (!lines.length) return null;
  const c0 = [T.R, 0], side = L => { const n = lineNormal(L); return Math.sign((c0[0] - L.q[0]) * n[0] + (c0[1] - L.q[1]) * n[1]) || 1; };
  const toTorus = (c, rho) => (c[0] > rho ? { ...T, c: add(T.c, mul(T.a, c[1])), R: c[0], r: rho } : null);
  if (lines.length === 1) {
    const L = lines[0], n = lineNormal(L), s = side(L);
    return rho => { const h = (c0[0] - L.q[0]) * n[0] + (c0[1] - L.q[1]) * n[1], k = s * rho - h; return toTorus([c0[0] + n[0] * k, c0[1] + n[1] * k], rho); };
  }
  // Dos rectas no paralelas: n1·(c - q1) = s1 ρ, n2·(c - q2) = s2 ρ.
  const [L1, L2] = lines, n1 = lineNormal(L1), n2 = lineNormal(L2), s1 = side(L1), s2 = side(L2);
  const det = n1[0] * n2[1] - n1[1] * n2[0]; if (Math.abs(det) < 1e-6) return null;
  return rho => {
    const b1 = s1 * rho + n1[0] * L1.q[0] + n1[1] * L1.q[1], b2 = s2 * rho + n2[0] * L2.q[0] + n2[1] * L2.q[1];
    return toTorus([(b1 * n2[1] - b2 * n1[1]) / det, (n1[0] * b2 - n2[0] * b1) / det], rho);
  };
}
