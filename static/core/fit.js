// Ajuste de primitivas analíticas (cilindro, cono, esfera, toro) a puntos de la malla.

export const sub = (a, b) => [a[0] - b[0], a[1] - b[1], a[2] - b[2]];
export const add = (a, b) => [a[0] + b[0], a[1] + b[1], a[2] + b[2]];
export const mul = (a, s) => [a[0] * s, a[1] * s, a[2] * s];
export const dot = (a, b) => a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
export const cross = (a, b) => [a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]];
export const norm = a => Math.hypot(a[0], a[1], a[2]);
export const unit = a => { const l = norm(a); return l > 0 ? [a[0] / l, a[1] / l, a[2] / l] : [0, 0, 1]; };
/** Vector unitario perpendicular a `a`. */
export function perp(a) {
  const b = Math.abs(a[0]) < 0.6 ? [1, 0, 0] : Math.abs(a[1]) < 0.6 ? [0, 1, 0] : [0, 0, 1];
  return unit(cross(a, b));
}

/** Resuelve A x = b (A densa n×n, se modifica). Devuelve null si es singular. */
export function solve(A, b) {
  const n = b.length;
  for (let i = 0; i < n; i++) {
    let p = i;
    for (let r = i + 1; r < n; r++) if (Math.abs(A[r][i]) > Math.abs(A[p][i])) p = r;
    if (Math.abs(A[p][i]) < 1e-300) return null;
    [A[i], A[p]] = [A[p], A[i]]; [b[i], b[p]] = [b[p], b[i]];
    for (let r = i + 1; r < n; r++) {
      const f = A[r][i] / A[i][i]; if (f === 0) continue;
      for (let c = i; c < n; c++) A[r][c] -= f * A[i][c];
      b[r] -= f * b[i];
    }
  }
  const x = new Array(n);
  for (let i = n - 1; i >= 0; i--) {
    let s = b[i];
    for (let c = i + 1; c < n; c++) s -= A[i][c] * x[c];
    x[i] = s / A[i][i];
  }
  return x.every(Number.isFinite) ? x : null;
}

/** Autovalores/autovectores de una matriz simétrica 3×3 (Jacobi). Orden ascendente. */
export function eigSym3(m) {
  const a = [[m[0][0], m[0][1], m[0][2]], [m[1][0], m[1][1], m[1][2]], [m[2][0], m[2][1], m[2][2]]];
  const v = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
  for (let sweep = 0; sweep < 60; sweep++) {
    const off = Math.abs(a[0][1]) + Math.abs(a[0][2]) + Math.abs(a[1][2]);
    if (off < 1e-300) break;
    for (let p = 0; p < 2; p++) for (let q = p + 1; q < 3; q++) {
      if (Math.abs(a[p][q]) < 1e-300) continue;
      const theta = (a[q][q] - a[p][p]) / (2 * a[p][q]);
      const t = Math.sign(theta || 1) / (Math.abs(theta) + Math.sqrt(theta * theta + 1));
      const c = 1 / Math.sqrt(t * t + 1), s = t * c;
      for (let k = 0; k < 3; k++) { const x = a[k][p], y = a[k][q]; a[k][p] = c * x - s * y; a[k][q] = s * x + c * y; }
      for (let k = 0; k < 3; k++) { const x = a[p][k], y = a[q][k]; a[p][k] = c * x - s * y; a[q][k] = s * x + c * y; }
      for (let k = 0; k < 3; k++) { const x = v[k][p], y = v[k][q]; v[k][p] = c * x - s * y; v[k][q] = s * x + c * y; }
    }
  }
  const order = [0, 1, 2].sort((i, j) => a[i][i] - a[j][j]);
  return { values: order.map(i => a[i][i]), vectors: order.map(i => [v[0][i], v[1][i], v[2][i]]) };
}

// ---------------------------------------------------------------- distancias

/** Distancia con signo de x a la primitiva (positiva hacia la normal natural). */
export function dist(S, x) {
  switch (S.type) {
    case 'plane': return dot(S.n, x) - S.d;
    case 'sphere': return Math.hypot(x[0] - S.c[0], x[1] - S.c[1], x[2] - S.c[2]) - S.r;
    case 'cylinder': {
      const v = sub(x, S.p), h = dot(v, S.a);
      return Math.sqrt(Math.max(0, dot(v, v) - h * h)) - S.r;
    }
    case 'cone': {
      const v = sub(x, S.apex), h = dot(v, S.a);
      const rho = Math.sqrt(Math.max(0, dot(v, v) - h * h));
      return rho * Math.cos(S.alpha) - h * Math.sin(S.alpha);
    }
    case 'torus': {
      const v = sub(x, S.c), h = dot(v, S.a);
      const rho = Math.sqrt(Math.max(0, dot(v, v) - h * h));
      return Math.hypot(rho - S.R, h) - S.r;
    }
  }
  throw new Error('tipo');
}

/** Normal natural de la primitiva en la proyección de x. */
export function normalAt(S, x) {
  switch (S.type) {
    case 'plane': return S.n;
    case 'sphere': return unit(sub(x, S.c));
    case 'cylinder': { const v = sub(x, S.p), h = dot(v, S.a); return unit(sub(v, mul(S.a, h))); }
    case 'cone': {
      const v = sub(x, S.apex), h = dot(v, S.a), rad = unit(sub(v, mul(S.a, h)));
      return unit(sub(mul(rad, Math.cos(S.alpha)), mul(S.a, Math.sin(S.alpha))));
    }
    case 'torus': {
      const v = sub(x, S.c), h = dot(v, S.a), radial = sub(v, mul(S.a, h)), rho = norm(radial), rad = unit(radial);
      return unit(add(mul(rad, rho - S.R), mul(S.a, h)));
    }
  }
  throw new Error('tipo');
}

/** Proyección de x sobre la primitiva. */
export function project(S, x) { return sub(x, mul(normalAt(S, x), dist(S, x))); }

// ------------------------------------------------------------ ajuste no lineal

const toVec = {
  sphere: S => [...S.c, S.r],
  cylinder: S => [...S.p, ...S.a, S.r],
  cone: S => [...S.apex, ...S.a, S.alpha],
  torus: S => [...S.c, ...S.a, S.R, S.r],
};
const fromVec = {
  sphere: v => ({ type: 'sphere', c: [v[0], v[1], v[2]], r: Math.abs(v[3]) }),
  cylinder: (v, ref) => {
    const a = unit([v[3], v[4], v[5]]); let p = [v[0], v[1], v[2]];
    const o = ref || p; p = add(p, mul(a, dot(sub(o, p), a)));
    return { type: 'cylinder', p, a, r: Math.abs(v[6]) };
  },
  cone: v => ({ type: 'cone', apex: [v[0], v[1], v[2]], a: unit([v[3], v[4], v[5]]), alpha: v[6] }),
  torus: v => ({ type: 'torus', c: [v[0], v[1], v[2]], a: unit([v[3], v[4], v[5]]), R: Math.abs(v[6]), r: Math.abs(v[7]) }),
};
const isDir = { sphere: [], cylinder: [3, 4, 5], cone: [3, 4, 5, 6], torus: [3, 4, 5] };

/** Refinamiento Levenberg–Marquardt de una primitiva sobre los puntos P. */
export function refine(S0, P, scale, iterations = 25) {
  const type = S0.type; if (type === 'plane') return S0;
  const centroid = [0, 0, 0];
  for (const p of P) { centroid[0] += p[0] / P.length; centroid[1] += p[1] / P.length; centroid[2] += p[2] / P.length; }
  const make = v => fromVec[type](v, centroid);
  const cost = S => { let s = 0; for (const p of P) { const d = dist(S, p); s += d * d; } return s; };
  let S = make(toVec[type](S0)), x = toVec[type](S), c = cost(S), lambda = 1e-3;
  const n = x.length, dirs = isDir[type];
  for (let it = 0; it < iterations; it++) {
    const r0 = P.map(p => dist(S, p));
    const J = [];
    for (let k = 0; k < n; k++) {
      const h = dirs.includes(k) ? 1e-6 : 1e-6 * scale;
      const xp = x.slice(); xp[k] += h; const xm = x.slice(); xm[k] -= h;
      const Sp = make(xp), Sm = make(xm);
      J.push(P.map(p => (dist(Sp, p) - dist(Sm, p)) / (2 * h)));
    }
    const JTJ = [], g = [];
    for (let i = 0; i < n; i++) {
      JTJ.push(new Array(n)); let gi = 0;
      for (let m = 0; m < P.length; m++) gi += J[i][m] * r0[m];
      g.push(-gi);
      for (let j = 0; j <= i; j++) { let s = 0; for (let m = 0; m < P.length; m++) s += J[i][m] * J[j][m]; JTJ[i][j] = s; }
    }
    for (let i = 0; i < n; i++) for (let j = i + 1; j < n; j++) JTJ[i][j] = JTJ[j][i];
    let improved = false;
    for (let tries = 0; tries < 8; tries++) {
      const M = JTJ.map((row, i) => row.map((v, j) => (i === j ? v * (1 + lambda) + 1e-12 * lambda : v)));
      const step = solve(M, g.slice());
      if (step) {
        const xn = x.map((v, i) => v + step[i]), Sn = make(xn), cn = cost(Sn);
        if (cn < c) {
          const rel = (c - cn) / (c || 1);
          S = Sn; x = toVec[type](Sn); c = cn; lambda = Math.max(lambda / 10, 1e-12); improved = true;
          if (rel < 1e-10) it = iterations;
          break;
        }
      }
      lambda *= 10;
    }
    if (!improved) break;
  }
  if (type === 'cone') S = normalizeCone(S, P);
  return S;
}

function normalizeCone(S, P) {
  // Eje orientado hacia la hoja que contiene los puntos; semiángulo en (0, π/2).
  let alpha = S.alpha % Math.PI; if (alpha < 0) alpha += Math.PI;
  let a = S.a;
  if (alpha > Math.PI / 2) { alpha = Math.PI - alpha; a = mul(a, -1); }
  let side = 0; for (const p of P) side += dot(sub(p, S.apex), a);
  if (side < 0) a = mul(a, -1);
  return { type: 'cone', apex: S.apex, a, alpha };
}

// ------------------------------------------------------------ ajustes iniciales

/** Esfera por mínimos cuadrados algebraicos. */
export function fitSphere(P) {
  if (P.length < 4) return null;
  const c0 = P[0];
  const A = [[0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], b = [0, 0, 0, 0];
  for (const q of P) {
    const p = sub(q, c0), row = [2 * p[0], 2 * p[1], 2 * p[2], 1], rhs = dot(p, p);
    for (let i = 0; i < 4; i++) { b[i] += row[i] * rhs; for (let j = 0; j < 4; j++) A[i][j] += row[i] * row[j]; }
  }
  const x = solve(A, b); if (!x) return null;
  const r2 = x[3] + x[0] * x[0] + x[1] * x[1] + x[2] * x[2];
  if (!(r2 > 0)) return null;
  return { type: 'sphere', c: add([x[0], x[1], x[2]], c0), r: Math.sqrt(r2) };
}

/** Círculo 2D por mínimos cuadrados algebraicos: devuelve [cx, cy, r]. */
export function fitCircle2D(pts) {
  if (pts.length < 3) return null;
  let mx = 0, my = 0; for (const p of pts) { mx += p[0] / pts.length; my += p[1] / pts.length; }
  const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], b = [0, 0, 0];
  for (const q of pts) {
    const x = q[0] - mx, y = q[1] - my, row = [2 * x, 2 * y, 1], rhs = x * x + y * y;
    for (let i = 0; i < 3; i++) { b[i] += row[i] * rhs; for (let j = 0; j < 3; j++) A[i][j] += row[i] * row[j]; }
  }
  const s = solve(A, b); if (!s) return null;
  const r2 = s[2] + s[0] * s[0] + s[1] * s[1];
  if (!(r2 > 0)) return null;
  return [s[0] + mx, s[1] + my, Math.sqrt(r2)];
}

/** Cilindro: eje a partir de las normales (todas perpendiculares al eje) y círculo proyectado. */
export function fitCylinder(P, Nrm) {
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const n of Nrm) for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) M[i][j] += n[i] * n[j];
  const e = eigSym3(M);
  if (e.values[1] < 1e-9 * e.values[2]) return null; // normales paralelas: es un plano
  return cylinderWithAxis(P, e.vectors[0]);
}
export function cylinderWithAxis(P, a) {
  const e1 = perp(a), e2 = cross(a, e1), o = P[0];
  const c = fitCircle2D(P.map(p => { const v = sub(p, o); return [dot(v, e1), dot(v, e2)]; }));
  if (!c) return null;
  return { type: 'cylinder', p: add(o, add(mul(e1, c[0]), mul(e2, c[1]))), a, r: c[2] };
}

/**
 * Eje de revolución por geometría de rectas: las normales de una superficie de
 * revolución cortan todas al eje. Devuelve { a, p } o null.
 */
export function revolutionAxis(C, Nrm) {
  const n = C.length; if (n < 6) return null;
  const o = [0, 0, 0]; for (const c of C) { o[0] += c[0] / n; o[1] += c[1] / n; o[2] += c[2] / n; }
  let scale = 0; for (const c of C) scale = Math.max(scale, norm(sub(c, o)));
  if (!(scale > 0)) return null;
  const A = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], B = [[0, 0, 0], [0, 0, 0], [0, 0, 0]], D = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let k = 0; k < n; k++) {
    const l = Nrm[k], m = cross(mul(sub(C[k], o), 1 / scale), l);
    for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) { A[i][j] += l[i] * l[j]; B[i][j] += l[i] * m[j]; D[i][j] += m[i] * m[j]; }
  }
  // Minimizar sobre cbar: cbar = -A^-1 B c ; queda (D - B^T A^-1 B) c.
  const X = [];
  for (let j = 0; j < 3; j++) {
    const col = solve(A.map(r => r.slice()), [B[0][j], B[1][j], B[2][j]]);
    if (!col) return null; X.push(col);
  }
  const Sc = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) {
    let s = D[i][j]; for (let k = 0; k < 3; k++) s -= B[k][i] * X[j][k];
    Sc[i][j] = s;
  }
  for (let i = 0; i < 3; i++) for (let j = i + 1; j < 3; j++) Sc[i][j] = Sc[j][i] = (Sc[i][j] + Sc[j][i]) / 2;
  const e = eigSym3(Sc), c = e.vectors[0];
  const cbar = [0, 0, 0];
  for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) cbar[i] -= X[j][i] * c[j];
  const x0 = cross(c, cbar);
  return { a: c, p: add(o, mul(x0, scale)) };
}

/** Coordenadas de perfil (altura h sobre el eje, radio rho). */
function profile(P, axis) {
  return P.map(p => { const v = sub(p, axis.p), h = dot(v, axis.a); return [h, Math.sqrt(Math.max(0, dot(v, v) - h * h))]; });
}

export function fitCone(P, C, Nrm) {
  const axis = revolutionAxis(C, Nrm); if (!axis) return null;
  const q = profile(P, axis);
  // Recta rho = m h + k en el plano de perfil.
  let sh = 0, sr = 0, shh = 0, shr = 0; const n = q.length;
  for (const [h, r] of q) { sh += h; sr += r; shh += h * h; shr += h * r; }
  const den = n * shh - sh * sh; if (Math.abs(den) < 1e-300) return null;
  const m = (n * shr - sh * sr) / den, k = (sr - m * sh) / n;
  if (Math.abs(m) < 1e-6) return null;
  const hApex = -k / m;
  let a = axis.a; if (m < 0) a = mul(a, -1);
  return { type: 'cone', apex: add(axis.p, mul(axis.a, hApex)), a, alpha: Math.atan(Math.abs(m)) };
}

export function fitTorus(P, C, Nrm) {
  const axis = revolutionAxis(C, Nrm); if (!axis) return null;
  const c = fitCircle2D(profile(P, axis)); if (!c) return null;
  if (c[1] <= 0) return null;
  return { type: 'torus', c: add(axis.p, mul(axis.a, c[0])), a: axis.a, R: c[1], r: c[2] };
}

/** Plano por mínimos cuadrados. Devuelve { n, d, o }. */
export function fitPlane(P) {
  const n = P.length, o = [0, 0, 0];
  for (const p of P) { o[0] += p[0] / n; o[1] += p[1] / n; o[2] += p[2] / n; }
  const M = [[0, 0, 0], [0, 0, 0], [0, 0, 0]];
  for (const p of P) { const v = sub(p, o); for (let i = 0; i < 3; i++) for (let j = 0; j < 3; j++) M[i][j] += v[i] * v[j]; }
  const e = eigSym3(M);
  return { type: 'plane', n: e.vectors[0], d: dot(e.vectors[0], o), o };
}

export const DOF = { plane: 3, cylinder: 5, cone: 6, sphere: 4, torus: 7 };
