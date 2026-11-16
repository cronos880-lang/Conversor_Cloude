// Ajuste de curvas a las cadenas de aristas que separan dos superficies.
import { dist, normalAt, fitPlane, fitCircle2D, dot, sub, add, mul, cross, unit, norm, perp, solve } from './fit.js';

/** Lleva x a la intersección de las superficies A y B (iteración de Newton). */
export function toIntersection(x, A, B, iterations = 12) {
  let p = x;
  for (let it = 0; it < iterations; it++) {
    const dA = dist(A, p), dB = dist(B, p);
    if (Math.abs(dA) + Math.abs(dB) < 1e-13) break;
    const nA = normalAt(A, p), nB = normalAt(B, p), c = dot(nA, nB), det = 1 - c * c;
    if (det < 1e-6) { p = sub(p, add(mul(nA, dA / 2), mul(nB, dB / 2))); continue; }
    const a = (-dA + c * dB) / det, b = (-dB + c * dA) / det;
    p = add(p, add(mul(nA, a), mul(nB, b)));
  }
  return p;
}

function lineFit(pts, tol) {
  const a = pts[0], b = pts[pts.length - 1], d = sub(b, a), len = norm(d);
  if (!(len > 0)) return null;
  const u = mul(d, 1 / len); let last = 0;
  for (let i = 1; i < pts.length - 1; i++) {
    const v = sub(pts[i], a), s = dot(v, u);
    if (norm(sub(v, mul(u, s))) > tol || s <= last || s >= len) return null;
    last = s;
  }
  return { type: 'line', p: a, d: u };
}

function circleFit(pts, closed, tol) {
  const n = closed ? pts.length - 1 : pts.length;
  if (n < 3) return null;
  const P = pts.slice(0, n), plane = fitPlane(P);
  const e1 = perp(plane.n), e2 = cross(plane.n, e1);
  for (const p of P) if (Math.abs(dot(plane.n, p) - plane.d) > tol) return null;
  const q = P.map(p => { const v = sub(p, plane.o); return [dot(v, e1), dot(v, e2)]; });
  const c = fitCircle2D(q); if (!c) return null;
  for (const [x, y] of q) if (Math.abs(Math.hypot(x - c[0], y - c[1]) - c[2]) > tol) return null;
  // Sentido de giro y monotonía angular.
  const ang = q.map(([x, y]) => Math.atan2(y - c[1], x - c[0]));
  let total = 0, sign = 0;
  const count = closed ? n : n - 1;
  for (let i = 0; i < count; i++) {
    let d = ang[(i + 1) % n] - ang[i];
    while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI;
    if (d === 0) return null;
    if (sign === 0) sign = Math.sign(d); else if (Math.sign(d) !== sign) return null;
    total += d;
  }
  if (closed && Math.abs(Math.abs(total) - 2 * Math.PI) > 1e-6) return null;
  if (!closed && Math.abs(total) >= 2 * Math.PI) return null;
  const axis = sign > 0 ? plane.n : mul(plane.n, -1);
  const center = add(plane.o, add(mul(e1, c[0]), mul(e2, c[1])));
  return { type: 'circle', c: center, a: axis, x: unit(sub(pts[0], add(center, mul(axis, dot(sub(pts[0], center), axis))))), r: c[2], sweep: Math.abs(total) };
}

/** Punto de una curva analítica en el parámetro angular t (para muestrear). */
export function circlePoint(C, t) {
  const y = cross(C.a, C.x);
  return add(C.c, add(mul(C.x, C.r * Math.cos(t)), mul(y, C.r * Math.sin(t))));
}

function tangents(pts, closed, A, B) {
  const n = pts.length, out = [];
  for (let i = 0; i < n; i++) {
    let prev = pts[i - 1], next = pts[i + 1];
    if (closed) { if (i === 0) prev = pts[n - 2]; if (i === n - 1) next = pts[1]; }
    let fd;
    if (prev && next) {
      // Tangente de la circunferencia que pasa por los tres puntos.
      const a = sub(pts[i], prev), b = sub(next, pts[i]), la = norm(a), lb = norm(b);
      fd = unit(add(mul(a, lb / la), mul(b, la / lb)));
    } else fd = null;
    let t = null;
    if (A && B) {
      const c = cross(normalAt(A, pts[i]), normalAt(B, pts[i]));
      if (norm(c) > 0.08) {
        t = unit(c);
        const ref = fd || unit(sub(next || pts[i], prev || pts[i]));
        if (dot(t, ref) < 0) t = mul(t, -1);
        if (dot(t, ref) < 0.5) t = null;
      }
    }
    out.push(t || fd);
  }
  // Extremos sin tangente: reflejar la del vecino respecto a la cuerda.
  for (const [i, j] of [[0, 1], [n - 1, n - 2]]) if (!out[i]) {
    const chord = unit(sub(pts[i === 0 ? 1 : n - 1], pts[i === 0 ? 0 : n - 2]));
    const tj = out[j] || chord;
    out[i] = unit(sub(mul(chord, 2 * dot(tj, chord)), tj));
  }
  return out;
}

function splineFit(pts, closed, A, B) {
  const T = tangents(pts, closed, A, B), ctrl = [pts[0]], knots = [0], mults = [4];
  let u = 0;
  for (let i = 0; i + 1 < pts.length; i++) {
    const chordV = sub(pts[i + 1], pts[i]), chord = norm(chordV), cu = mul(chordV, 1 / chord);
    const h0 = (2 / 3) * chord / (1 + Math.max(0.2, dot(T[i], cu))), h1 = (2 / 3) * chord / (1 + Math.max(0.2, dot(T[i + 1], cu)));
    ctrl.push(add(pts[i], mul(T[i], h0)), sub(pts[i + 1], mul(T[i + 1], h1)), pts[i + 1]);
    u += chord; knots.push(u); mults.push(3);
  }
  mults[mults.length - 1] = 4;
  return { type: 'bspline', degree: 3, ctrl, knots, mults };
}

/** Punto de la B-spline cúbica (formada por tramos de Bézier) a mitad del tramo i. */
function splineMid(C, i) {
  const [a, b, c, d] = C.ctrl.slice(i * 3, i * 3 + 4);
  return [0, 1, 2].map(k => (a[k] + 3 * b[k] + 3 * c[k] + d[k]) / 8);
}

/**
 * Curva para una cadena de puntos entre las superficies A y B.
 * Devuelve la curva y `dev`, la desviación máxima medida respecto a ambas superficies.
 */
export function fitCurve(pts, closed, A, B, tol) {
  const onBoth = p => Math.max(Math.abs(dist(A, p)), Math.abs(dist(B, p)));
  let curve = null, dev = 0;
  if (!closed) {
    curve = lineFit(pts, tol);
    if (curve) {
      for (let i = 0; i + 1 < pts.length; i++) dev = Math.max(dev, onBoth(mul(add(pts[i], pts[i + 1]), 0.5)));
      if (dev <= 4 * tol || pts.length === 2) return { curve, dev };
      curve = null;
    }
  }
  curve = circleFit(pts, closed, tol);
  if (curve) {
    dev = 0;
    const y = cross(curve.a, curve.x), n = pts.length;
    for (let i = 0; i + 1 < n; i++) {
      const m = sub(mul(add(pts[i], pts[i + 1]), 0.5), curve.c);
      const dir = unit(sub(m, mul(curve.a, dot(m, curve.a))));
      // Punto medio del arco entre dos vértices consecutivos (el arco menor de cada tramo).
      dev = Math.max(dev, onBoth(add(curve.c, mul(dir, curve.r))));
      void y;
    }
    if (dev <= 4 * tol) return { curve, dev };
  }
  curve = splineFit(pts, closed, A, B);
  dev = 0;
  for (let i = 0; i + 1 < pts.length; i++) dev = Math.max(dev, onBoth(splineMid(curve, i)));
  return { curve, dev };
}
