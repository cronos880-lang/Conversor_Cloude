// Ajuste de curvas a las cadenas de aristas que separan dos superficies.
import { coaxialCircles } from './revolve.js';
import { dist, normalAt, project, fitPlane, fitCircle2D, dot, sub, add, mul, cross, unit, norm, perp, solve } from './fit.js';

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
  let nrm = plane.n;
  const first = pts[0], last = pts[pts.length - 1];
  let chord = null;
  if (!closed) {
    // El plano de la circunferencia debe contener la cuerda entre los extremos.
    const d = sub(last, first), len = norm(d); if (!(len > 0)) return null;
    chord = mul(d, 1 / len);
    nrm = sub(nrm, mul(chord, dot(nrm, chord)));
    if (norm(nrm) < 0.5) return null;
    nrm = unit(nrm);
  }
  const o = closed ? plane.o : first;
  const e1 = closed ? perp(nrm) : chord, e2 = cross(nrm, e1);
  for (const p of P) if (Math.abs(dot(nrm, sub(p, o))) > tol) return null;
  const q = P.map(p => { const v = sub(p, o); return [dot(v, e1), dot(v, e2)]; });
  let c = fitCircle2D(q); if (!c) return null;
  if (closed) c = [c[0], c[1], Math.hypot(q[0][0] - c[0], q[0][1] - c[1])];
  else { const half = q[n - 1][0] / 2; c = [half, c[1], Math.hypot(half, c[1])]; }
  for (const [x, y] of q) if (Math.abs(Math.hypot(x - c[0], y - c[1]) - c[2]) > tol) return null;
  // Sentido de giro y monotonía angular.
  const ang = q.map(([x, y]) => Math.atan2(y - c[1], x - c[0]));
  let total = 0, sign = 0, widest = 0;
  const count = closed ? n : n - 1;
  for (let i = 0; i < count; i++) {
    let d = ang[(i + 1) % n] - ang[i];
    while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI;
    if (d === 0) return null;
    widest = Math.max(widest, Math.abs(d));
    if (sign === 0) sign = Math.sign(d); else if (Math.sign(d) !== sign) return null;
    total += d;
  }
  // Un arco apoyado en dos grupos de puntos separados por un gran vacío no es fiable.
  if (count >= 3 && widest > 0.6 * Math.abs(total)) return null;
  if (closed && Math.abs(Math.abs(total) - 2 * Math.PI) > 1e-6) return null;
  if (!closed && Math.abs(total) >= 2 * Math.PI - 1e-9) return null;
  const axis = sign > 0 ? nrm : mul(nrm, -1);
  const center = add(o, add(mul(e1, c[0]), mul(e2, c[1])));
  return { type: 'circle', c: center, a: axis, x: unit(sub(first, center)), r: c[2], sweep: Math.abs(total) };
}


/**
 * Circunferencia exacta de la intersección entre una superficie de revolución y un plano
 * perpendicular a su eje, o entre dos superficies de revolución coaxiales tangentes.
 * Devuelve candidatos { c, a, r } (puede haber dos en el toro).
 */
function analyticCircles(A, B) {
  const known = analyticCirclesBasic(A, B);
  if (known.length || globalThis.process?.env?.NOCOAX) return known;
  // Cualquier otra pareja coaxial (cono–toro, cono–cilindro, toro–toro...): corte de perfiles.
  return coaxialCircles(A, B);
}

function analyticCirclesBasic(A, B) {
  const rev = S => S.type === 'cylinder' || S.type === 'cone' || S.type === 'torus' || S.type === 'sphere';
  const axisPt = S => S.p || S.apex || S.c;
  if (B.type === 'plane' && rev(A)) [A, B] = [B, A];
  const out = [];
  if (A.type === 'plane' && rev(B)) {
    const n = A.n, a = B.type === 'sphere' ? n : B.a;
    if (Math.abs(dot(n, a)) < 1 - 1e-9) return out;
    const o = axisPt(B), t = (A.d - dot(n, o)) / dot(n, a), c = add(o, mul(a, t));
    if (B.type === 'cylinder') out.push(B.r);
    else if (B.type === 'cone') { const h = dot(sub(c, B.apex), B.a); if (h > 0) out.push(h * Math.tan(B.alpha)); }
    else if (B.type === 'sphere') { const h = dot(n, B.c) - A.d, q = B.r * B.r - h * h; out.push(Math.sqrt(Math.max(0, q))); }
    else {
      const h = dot(n, B.c) - A.d, q = B.r * B.r - h * h, w = Math.sqrt(Math.max(0, q));
      // Plano tangente al tubo: el borde puede quedar en cualquier punto de la franja de contacto.
      if (w < 1e-6) return [{ c, a, r: B.R, tangent: B.r }];
      out.push(B.R - w, B.R + w);
    }
    return out.filter(r => r > 0).map(r => ({ c, a, r }));
  }
  if (rev(A) && rev(B) && A.type !== 'sphere' && B.type !== 'sphere') {
    if (Math.abs(dot(A.a, B.a)) < 1 - 1e-9) return out;
    const T = A.type === 'torus' ? A : B.type === 'torus' ? B : null, O = T === A ? B : A;
    if (!T || O.type !== 'cylinder') return out;
    const v = sub(T.c, O.p); if (norm(sub(v, mul(O.a, dot(v, O.a)))) > 1e-6) return out;
    // Toro tangente al cilindro coaxial: la circunferencia de contacto está en el plano del centro.
    if (Math.abs(Math.abs(T.R - O.r) - T.r) < 1e-6) out.push({ c: T.c, a: T.a, r: O.r, tangent: T.r });
  }
  return out;
}

/** Rectas exactas de la intersección de un plano con un cilindro de eje paralelo al plano. */
function analyticLines(A, B, tol = 0) {
  if (B.type === 'plane' && A.type === 'cylinder') [A, B] = [B, A];
  if (A.type !== 'plane' || B.type !== 'cylinder' || Math.abs(dot(A.n, B.a)) > 1e-9) return [];
  const h = dot(A.n, B.p) - A.d, q = B.r * B.r - h * h; if (q < 0) return [];
  const foot = sub(B.p, mul(A.n, h)), side = cross(A.n, B.a), w = Math.sqrt(q);
  // Corte rasante: el borde de la malla puede caer en una franja ancha alrededor de la recta.
  const band = Math.max(3 * tol, Math.min(2 * tol * B.r / Math.max(w, 1e-12), 2 * Math.sqrt(2 * B.r * tol)));
  return [{ p: add(foot, mul(side, w)), d: B.a, band }, { p: sub(foot, mul(side, w)), d: B.a, band }];
}

/** Arco o circunferencia { c, a, r } que pasa por los puntos (dentro de lim), o null. */
function circleOn(C, pts, closed, lim, back = 0) {
  const x0 = sub(pts[0], C.c), xr = sub(x0, mul(C.a, dot(x0, C.a))); if (norm(xr) < 1e-9) return null;
  const ex = unit(xr), ey = cross(C.a, ex);
  const n = closed ? pts.length - 1 : pts.length, ang = [];
  for (let i = 0; i < n; i++) {
    const v = sub(pts[i], C.c), h = dot(v, C.a), rad = Math.hypot(dot(v, ex), dot(v, ey));
    if (Math.hypot(h, rad - C.r) > lim) return null;
    ang.push(Math.atan2(dot(v, ey), dot(v, ex)));
  }
  let total = 0, sign = 0, widest = 0; const count = closed ? n : n - 1;
  for (let i = 0; i < count; i++) {
    let d = ang[(i + 1) % n] - ang[i];
    while (d > Math.PI) d -= 2 * Math.PI; while (d < -Math.PI) d += 2 * Math.PI;
    if (d === 0) return null;
    widest = Math.max(widest, Math.abs(d));
    // En contactos tangentes la cadena de la malla puede retroceder un poco en ángulo.
    if (Math.abs(d) <= back) { total += d; continue; }
    if (sign === 0) sign = Math.sign(d); else if (Math.sign(d) !== sign) return null;
    total += d;
  }
  if (sign === 0) sign = Math.sign(total) || 1;
  if (sign * total <= 0) return null;
  if (count >= 3 && widest > 0.6 * Math.abs(total)) return null;
  if (closed && Math.abs(Math.abs(total) - 2 * Math.PI) > 1e-6) return null;
  if (!closed && Math.abs(total) >= 2 * Math.PI - 1e-9) return null;
  const axis = sign > 0 ? C.a : mul(C.a, -1);
  return { type: 'circle', c: C.c, a: axis, x: ex, r: C.r, sweep: Math.abs(total) };
}

/** Punto común a A y B también cuando son tangentes (proyecciones alternadas). */
function onIntersection(x, A, B) {
  let q = toIntersection(x, A, B);
  const err = p => Math.max(Math.abs(dist(A, p)), Math.abs(dist(B, p)));
  if (q.every(Number.isFinite) && err(q) < 1e-7) return q;
  q = x;
  for (let it = 0; it < 40 && err(q) > 1e-7; it++) q = project(B, project(A, q));
  return q;
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

/** Lleva p al punto más cercano común a varias superficies (mínimos cuadrados). */
export function snapPoint(p, surfaces, maxMove) {
  let x = p;
  for (let it = 0; it < 12; it++) {
    const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], g = [0, 0, 0]; let err = 0;
    for (const S of surfaces) {
      const d = dist(S, x), n = normalAt(S, x); err = Math.max(err, Math.abs(d));
      for (let i = 0; i < 3; i++) { g[i] -= n[i] * d; for (let j = 0; j < 3; j++) M[i][j] += n[i] * n[j]; }
    }
    if (err < 1e-13) break;
    for (let i = 0; i < 3; i++) M[i][i] += 1e-9 * surfaces.length;
    const step = solve(M, g); if (!step) break;
    x = add(x, step);
  }
  if (!x.every(Number.isFinite) || norm(sub(x, p)) > maxMove) return p;
  return x;
}

/**
 * Curva para una cadena de puntos entre las superficies A y B. Los extremos se
 * respetan exactamente. Devuelve la curva, `dev` (desviación máxima medida
 * respecto a ambas superficies) y `ok` (la curva sigue de cerca a la cadena).
 */
export function fitCurve(pts, closed, A, B, tol) {
  const onBoth = p => Math.max(Math.abs(dist(A, p)), Math.abs(dist(B, p)));
  // La curva no puede alejarse de la poligonal más que la flecha de un arco razonable.
  const near = (mid, i) => { const c = sub(pts[i + 1], pts[i]); return norm(sub(mid, mul(add(pts[i], pts[i + 1]), 0.5))) <= Math.max(3 * tol, 0.13 * norm(c)); };
  let curve = null, dev = 0;
  if (!closed && !globalThis.process?.env?.NOANALYTIC) {
    for (const L of analyticLines(A, B, tol)) {
      const off = p => { const v = sub(p, L.p); return norm(sub(v, mul(L.d, dot(v, L.d)))); };
      if (off(pts[0]) > 0.02 * tol || off(pts[pts.length - 1]) > 0.02 * tol || pts.some(p => off(p) > L.band)) continue;
      const d = sub(pts[pts.length - 1], pts[0]), len = norm(d); if (!(len > 0)) continue;
      return { curve: { type: 'line', p: pts[0], d: mul(d, 1 / len) }, dev: 0, ok: true };
    }
  }
  if (!closed) {
    curve = lineFit(pts, tol);
    if (curve) {
      for (let i = 0; i + 1 < pts.length; i++) dev = Math.max(dev, onBoth(mul(add(pts[i], pts[i + 1]), 0.5)));
      if (dev <= 4 * tol) return { curve, dev, ok: true };
      if (pts.length === 2) {
        // Tramo de dos puntos sobre una superficie curva: puntos intermedios de la intersección real.
        const dense = [pts[0]];
        for (let k = 1; k < 8; k++) {
          const q = onIntersection(add(mul(pts[0], 1 - k / 8), mul(pts[1], k / 8)), A, B);
          if (!(onBoth(q) < 1e-6) || norm(sub(q, pts[0])) > 2 * norm(sub(pts[1], pts[0]))) return { curve, dev, ok: true };
          dense.push(q);
        }
        dense.push(pts[1]);
        return fitCurve(dense, false, A, B, tol);
      }
      curve = null;
    }
  }
  // Intersección conocida analíticamente (redondeos tangentes, taladros perpendiculares...).
  for (const C of globalThis.process?.env?.NOANALYTIC ? [] : analyticCircles(A, B)) {
    // Los extremos (vértices ya calculados) deben estar sobre la circunferencia exacta.
    const offC = p => { const v = sub(p, C.c), h = dot(v, C.a); return Math.hypot(h, norm(sub(v, mul(C.a, h))) - C.r); };
    if (offC(pts[0]) > 0.02 * tol || offC(pts[pts.length - 1]) > 0.02 * tol) continue;
    const c = circleOn(C, pts, closed, C.tangent ? Math.max(3 * tol, 2 * Math.sqrt(2 * C.tangent * tol)) : 3 * tol, C.tangent ? 0.08 : 0);
    if (c && (closed || pts.length > 2 || c.sweep > 0.02)) return { curve: c, dev: 0, ok: true };
  }
  curve = circleFit(pts, closed, tol);
  if (curve) {
    dev = 0; let ok = true;
    for (let i = 0; i + 1 < pts.length; i++) {
      const m = sub(mul(add(pts[i], pts[i + 1]), 0.5), curve.c);
      const mid = add(curve.c, mul(unit(sub(m, mul(curve.a, dot(m, curve.a)))), curve.r));
      dev = Math.max(dev, onBoth(mid)); ok = ok && near(mid, i);
    }
    if (dev <= 4 * tol && ok) return { curve, dev, ok };
  }
  // Spline por los puntos; si se separa de las superficies, densificar con puntos de la intersección.
  let ok = true;
  { const first = splineFit(pts, closed, A, B); for (let i = 0; i + 1 < pts.length; i++) ok = ok && near(splineMid(first, i), i); }
  let P = pts;
  for (let round = 0; ; round++) {
    curve = splineFit(P, closed, A, B);
    dev = 0; const mids = [];
    for (let i = 0; i + 1 < P.length; i++) { const m = splineMid(curve, i); mids.push(m); dev = Math.max(dev, onBoth(m)); }
    if (!ok || dev < 1e-5 || round === 6 || P.length > 3000) break;
    const next = [P[0]]; let moved = false;
    for (let i = 0; i + 1 < P.length; i++) {
      if (onBoth(mids[i]) > 1e-5) {
        const q = onIntersection(mids[i], A, B);
        if (onBoth(q) < 1e-6 && norm(sub(q, mids[i])) < 0.3 * norm(sub(P[i + 1], P[i]))) { next.push(q); moved = true; }
      }
      next.push(P[i + 1]);
    }
    if (!moved) break;
    P = next;
  }
  return { curve, dev, ok };
}

/**
 * Curva suave de respaldo para cadenas que no siguen bien la intersección exacta: una sola
 * B-spline por puntos espaciados al menos `gap`, sin vértices intermedios (que crearían
 * aristas diminutas). Solo se acepta si se mantiene cerca de ambas superficies.
 */
export function smoothCurve(pts, closed, A, B, tol, gap) {
  if (pts.length < 3) return null;
  const last = pts[pts.length - 1], kept = [pts[0]];
  for (let i = 1; i < pts.length - 1; i++) if (norm(sub(pts[i], kept[kept.length - 1])) >= gap && norm(sub(pts[i], last)) >= gap) kept.push(pts[i]);
  kept.push(last);
  if (kept.length < 2 || norm(sub(kept[0], kept[1])) < 1e-9) return null;
  const curve = splineFit(kept, closed, null, null);
  const onBoth = p => Math.max(Math.abs(dist(A, p)), Math.abs(dist(B, p)));
  let dev = 0;
  for (let i = 0; i + 1 < kept.length; i++) {
    const m = splineMid(curve, i), c = mul(add(kept[i], kept[i + 1]), 0.5);
    if (norm(sub(m, c)) > Math.max(3 * tol, 0.13 * norm(sub(kept[i + 1], kept[i])))) return null;
    dev = Math.max(dev, onBoth(m));
  }
  for (const p of pts) dev = Math.max(dev, onBoth(p));
  return dev <= 4 * tol ? { curve, dev, ok: true } : null;
}

/**
 * Vértice sobre el contacto de dos superficies tangentes (p. ej. redondeo y plano): la
 * posición a lo ancho de la franja de contacto no queda definida por las superficies, así
 * que se toma el punto de la circunferencia de contacto que está sobre las demás caras.
 */
export function cornerOnTangent(p, surfaces, tol) {
  for (let i = 0; i < surfaces.length; i++) for (let j = i + 1; j < surfaces.length; j++) {
    for (const C of analyticCircles(surfaces[i], surfaces[j])) {
      if (!C.tangent) continue;
      const band = Math.max(3 * tol, 2 * Math.sqrt(2 * C.tangent * tol));
      const v = sub(p, C.c), h = dot(v, C.a), rv = sub(v, mul(C.a, h)), rl = norm(rv);
      if (rl < 1e-9 || Math.hypot(h, rl - C.r) > band) continue;
      const ex = mul(rv, 1 / rl), ey = cross(C.a, ex), at = t => add(C.c, add(mul(ex, C.r * Math.cos(t)), mul(ey, C.r * Math.sin(t))));
      const others = surfaces.filter((_, k) => k !== i && k !== j);
      let t = 0;
      if (others.length) {
        const f = t => others.reduce((s, S) => s + dist(S, at(t)) ** 2, 0);
        for (let it = 0; it < 30; it++) {
          const e = 1e-7, g = (f(t + e) - f(t - e)) / (2 * e), H = (f(t + e) - 2 * f(t) + f(t - e)) / (e * e);
          if (!(H > 0)) break;
          const step = g / H; t -= Math.max(-0.05, Math.min(0.05, step));
          if (Math.abs(step) < 1e-12) break;
        }
        if (!others.every(S => Math.abs(dist(S, at(t))) < Math.max(1e-6, 0.1 * tol))) continue;
      }
      const q = at(t);
      if (norm(sub(q, p)) <= band) return q;
    }
  }
  return null;
}
