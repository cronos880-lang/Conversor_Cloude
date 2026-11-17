// Topología B-rep a partir de la segmentación: caras, bucles, aristas y vértices.
import { fitCurve, snapPoint, toIntersection } from './curves.js';
import { dist, normalAt, cross as crossV, norm as normV } from './fit.js';
import { perp, cross, add, mul } from './fit.js';

/**
 * Devuelve { vertices, edges, faces, shells }.
 *  - vertices[i] = [x, y, z]
 *  - edges[i] = { v0, v1, curve, dev }
 *  - faces[i] = { S, sign, loops: [[{ edge, forward }]] , shell }
 *  - shells[i] = { faces: [índices], volume }
 */
export function buildBrep(mesh, seg, comp, options = {}) {
  const polyFaces = options.polyFaces || new Set();
  const { V, T, adj, nt } = mesh, { reg, prims, tol } = seg;
  const P = i => [V[i * 3], V[i * 3 + 1], V[i * 3 + 2]];
  const next = h => (h % 3 === 2 ? h - 2 : h + 1);
  const regOf = h => reg[(h / 3) | 0];
  const boundary = h => regOf(h) !== regOf(adj[h]);
  for (let h = 0; h < nt * 3; h++) if (adj[h] < 0) throw new Error('La malla no es cerrada.');

  // Vértices esquina: aquellos con un número de aristas frontera distinto de 2.
  const count = new Int32Array(mesh.nv);
  for (let h = 0; h < nt * 3; h++) if (boundary(h)) count[T[h]]++;

  const visited = new Uint8Array(nt * 3);
  const faces = prims.map((p, i) => ({ S: p.S, sign: p.sign, loops: [], shell: -1, region: i }));
  // 1) Topología: cadenas de aristas de malla entre esquinas.
  const chains = [], chainOf = new Map(); // primera semiarista -> { edge, forward }
  const cornerRegions = new Map();
  const compOfRegion = new Int32Array(prims.length);
  for (let t = 0; t < nt; t++) compOfRegion[reg[t]] = comp.comp[t];
  const touch = (v0, r) => { const v = v0 + compOfRegion[r] * mesh.nv; let s = cornerRegions.get(v); if (!s) cornerRegions.set(v, s = new Set()); s.add(r); };
  for (let h0 = 0; h0 < nt * 3; h0++) {
    if (visited[h0] || !boundary(h0)) continue;
    const A = regOf(h0), loop = [];
    let h = h0;
    do {
      visited[h] = 1; loop.push(h);
      let n = next(h);
      while (regOf(adj[n]) === A) n = next(adj[n]);
      h = n;
    } while (h !== h0);
    // Girar el bucle para empezar en una esquina (o en el vértice de menor índice).
    let start = -1;
    for (let i = 0; i < loop.length; i++) if (count[T[loop[i]]] !== 2) { start = i; break; }
    if (start < 0) { start = 0; for (let i = 1; i < loop.length; i++) if (T[loop[i]] < T[loop[start]]) start = i; }
    const rot = loop.slice(start).concat(loop.slice(0, start));
    const refs = [];
    for (let i = 0; i < rot.length;) {
      let j = i + 1;
      while (j < rot.length && count[T[rot[j]]] === 2) j++;
      const chain = rot.slice(i, j), first = chain[0];
      let ref = chainOf.get(first);
      if (!ref) {
        const B = regOf(adj[first]);
        const vs = chain.map(x => T[x]); vs.push(T[next(chain[chain.length - 1])]);
        const id = chains.length; chains.push({ vs, A, B });
        touch(vs[0], A); touch(vs[0], B); touch(vs[vs.length - 1], A); touch(vs[vs.length - 1], B);
        ref = { edge: id, forward: true };
        chainOf.set(adj[chain[chain.length - 1]], { edge: id, forward: false });
      }
      refs.push(ref); i = j;
    }
    faces[A].loops.push(refs);
  }
  // 2) Geometría: vértices y puntos intermedios llevados sobre las superficies ajustadas.
  const maxMove = 4 * tol;
  const vertices = [], vertexId = new Map();
  const vertexOf = (v0, shell) => {
    const v = v0 + shell * mesh.nv; let id = vertexId.get(v);
    if (id === undefined) {
      id = vertices.length; vertexId.set(v, id);
      // Junto a una faceta el vértice de la malla ya es exacto: ajustarlo a planos casi
      // paralelos lo desplazaría y podría cruzar caras.
      const regs = [...cornerRegions.get(v)];
      vertices.push(regs.some(r => prims[r].facet) ? P(v0) : snapPoint(P(v0), regs.map(r => prims[r].S), maxMove));
    }
    return id;
  };
  const edges = [], chainInfo = [];
  const same = (p, q) => Math.hypot(p[0] - q[0], p[1] - q[1], p[2] - q[2]) < 1e-9;
  const chainEdges = chains.map(({ vs, A, B }) => {
    const SA = prims[A].S, SB = prims[B].S, closed = vs[0] === vs[vs.length - 1];
    const v0 = vertexOf(vs[0], compOfRegion[A]), v1 = vertexOf(vs[vs.length - 1], compOfRegion[A]);
    let pts = vs.map((v, i) => {
      if (i === 0) return vertices[v0];
      if (i === vs.length - 1) return vertices[v1];
      const p = P(v);
      // Intersección transversal: el punto exacto común; tangente: punto intermedio.
      const transversal = normV(crossV(normalAt(SA, p), normalAt(SB, p))) > 0.12;
      const q = transversal ? toIntersection(p, SA, SB) : snapPoint(p, [SA, SB], maxMove);
      return q.every(Number.isFinite) && Math.hypot(q[0] - p[0], q[1] - p[1], q[2] - p[2]) <= maxMove ? q : p;
    });
    // Quitar puntos repetidos tras el ajuste.
    pts = pts.filter((p, i) => i === 0 || i === pts.length - 1 || (!same(p, pts[i - 1]) && !same(p, pts[pts.length - 1])));
    let fit = null;
    if (!polyFaces.has(A) && !polyFaces.has(B)) { fit = fitCurve(pts, closed, SA, SB, tol); if (!fit.ok) fit = null; }
    chainInfo.push({ A: SA.type, B: SB.type, kind: fit ? fit.curve.type : 'poly', n: pts.length, closed, len: +Math.hypot(...pts[0].map((x, i) => x - pts[pts.length - 1][i])).toFixed(3), resA: +(prims[A].res || 0).toExponential(1), resB: +(prims[B].res || 0).toExponential(1) });
    if (fit) { edges.push({ v0, v1, curve: fit.curve, dev: fit.dev, n: pts.length, A, B }); return [edges.length - 1]; }
    // Poligonal: un segmento recto por cada arista de la malla.
    const ids = []; let prev = v0;
    for (let i = 1; i < pts.length; i++) {
      let cur = v1;
      if (i < pts.length - 1) { cur = vertices.length; vertices.push(pts[i]); }
      const a = vertices[prev], b = vertices[cur], d = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], len = Math.hypot(d[0], d[1], d[2]);
      const mid = [(a[0] + b[0]) / 2, (a[1] + b[1]) / 2, (a[2] + b[2]) / 2];
      edges.push({ v0: prev, v1: cur, curve: { type: 'line', p: a, d: len > 0 ? [d[0] / len, d[1] / len, d[2] / len] : [0, 0, 1] }, dev: Math.max(Math.abs(dist(SA, mid)), Math.abs(dist(SB, mid))), n: 2, poly: true });
      ids.push(edges.length - 1); prev = cur;
    }
    return ids;
  });
  // En las caras planas, el contorno exterior es el bucle que encierra más área.
  for (const f of faces) {
    if (f.S.type !== 'plane' || f.loops.length < 2) continue;
    const n = f.S.n; let best = -Infinity;
    f.loops.forEach((loop, i) => {
      const pts = loop.flatMap(ref => { const vs = chains[ref.edge].vs.map(P); return (ref.forward ? vs : vs.slice().reverse()).slice(0, -1); });
      let ax = 0, ay = 0, az = 0;
      for (let k = 0; k < pts.length; k++) { const p = pts[k], q = pts[(k + 1) % pts.length]; ax += p[1] * q[2] - p[2] * q[1]; ay += p[2] * q[0] - p[0] * q[2]; az += p[0] * q[1] - p[1] * q[0]; }
      const area = (ax * n[0] + ay * n[1] + az * n[2]) / 2;
      if (area > best) { best = area; f.outer = i; }
    });
  }
  for (const f of faces) f.loops = f.loops.map(loop => loop.flatMap(ref => {
    const list = chainEdges[ref.edge];
    return ref.forward ? list.map(e => ({ edge: e, forward: true })) : list.slice().reverse().map(e => ({ edge: e, forward: false }));
  }));

  // Envolventes: una por componente conexa de la malla.
  const shells = [];
  for (let c = 0; c < comp.count; c++) shells.push({ faces: [], volume: 0 });
  const o = mesh.min;
  for (let t = 0; t < nt; t++) {
    const f = faces[reg[t]], c = comp.comp[t];
    if (f.shell < 0) { f.shell = c; shells[c].faces.push(reg[t]); }
    const a = P(T[t * 3]), b = P(T[t * 3 + 1]), d = P(T[t * 3 + 2]);
    const x = [a, b, d].map(p => [p[0] - o[0], p[1] - o[1], p[2] - o[2]]);
    shells[c].volume += (x[0][0] * (x[1][1] * x[2][2] - x[1][2] * x[2][1]) + x[0][1] * (x[1][2] * x[2][0] - x[1][0] * x[2][2]) + x[0][2] * (x[1][0] * x[2][1] - x[1][1] * x[2][0])) / 6;
  }

  // Cavidades: una envolvente de volumen negativo es un hueco del sólido que la contiene.
  for (let c = 0; c < shells.length; c++) {
    if (shells[c].volume >= 0) continue;
    let t0 = 0; while (comp.comp[t0] !== c) t0++;
    const p = P(T[t0 * 3]); let parent = -1;
    for (let k = 0; k < shells.length; k++) {
      if (shells[k].volume <= 0 || (parent >= 0 && shells[k].volume >= shells[parent].volume)) continue;
      // Rayo en +x (con una ligera inclinación para evitar aristas): número impar de cruces = dentro.
      const dir = [1, 0.0137, 0.0071]; let hits = 0;
      for (let t = 0; t < nt; t++) {
        if (comp.comp[t] !== k) continue;
        const a = P(T[t * 3]), e1 = P(T[t * 3 + 1]).map((x, i) => x - a[i]), e2 = P(T[t * 3 + 2]).map((x, i) => x - a[i]);
        const h = [dir[1] * e2[2] - dir[2] * e2[1], dir[2] * e2[0] - dir[0] * e2[2], dir[0] * e2[1] - dir[1] * e2[0]];
        const det = e1[0] * h[0] + e1[1] * h[1] + e1[2] * h[2]; if (Math.abs(det) < 1e-300) continue;
        const s = p.map((x, i) => x - a[i]), u = (s[0] * h[0] + s[1] * h[1] + s[2] * h[2]) / det; if (u < 0 || u > 1) continue;
        const q = [s[1] * e1[2] - s[2] * e1[1], s[2] * e1[0] - s[0] * e1[2], s[0] * e1[1] - s[1] * e1[0]];
        const v = (dir[0] * q[0] + dir[1] * q[1] + dir[2] * q[2]) / det; if (v < 0 || u + v > 1) continue;
        if ((e2[0] * q[0] + e2[1] * q[1] + e2[2] * q[2]) / det > 0) hits++;
      }
      if (hits % 2 === 1) parent = k;
    }
    shells[c].parent = parent;
  }

  // Superficies cerradas sin contorno (esfera o toro completos): dividir en dos mitades.
  const extra = [];
  for (const f of faces) {
    if (f.loops.length || f.shell < 0) continue;
    const S = f.S, up = f.sign > 0;
    if (S.type === 'sphere') {
      const a = [0, 0, 1], x = [1, 0, 0], v = vertices.length; vertices.push(add(S.c, mul(x, S.r)));
      const e = edges.length; edges.push({ v0: v, v1: v, curve: { type: 'circle', c: S.c, a, x, r: S.r }, dev: 0, n: 0 });
      f.loops.push([{ edge: e, forward: up }]);
      const g = { S, sign: f.sign, loops: [[{ edge: e, forward: !up }]], shell: f.shell, region: f.region };
      shells[f.shell].faces.push(faces.length + extra.length); extra.push(g);
    } else if (S.type === 'torus') {
      const x = perp(S.a), ids = [];
      for (const r of [S.R + S.r, S.R - S.r]) {
        const v = vertices.length; vertices.push(add(S.c, mul(x, r)));
        ids.push(edges.length); edges.push({ v0: v, v1: v, curve: { type: 'circle', c: S.c, a: S.a, x, r }, dev: 0, n: 0 });
      }
      f.loops.push([{ edge: ids[0], forward: up }], [{ edge: ids[1], forward: !up }]);
      const g = { S, sign: f.sign, loops: [[{ edge: ids[0], forward: !up }], [{ edge: ids[1], forward: up }]], shell: f.shell, region: f.region };
      shells[f.shell].faces.push(faces.length + extra.length); extra.push(g);
    } else throw new Error('Superficie sin contorno no admitida.' + (globalThis.DEBUG_BREP ? ' ' + JSON.stringify({ type: S.type, region: f.region, shell: f.shell, tris: Array.from(seg.reg).filter(x => x === f.region).length }) : ''));
  }
  void cross;
  return { vertices, edges, faces: faces.concat(extra), shells, chainInfo };
}
