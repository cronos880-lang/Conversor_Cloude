// Topología B-rep a partir de la segmentación: caras, bucles, aristas y vértices.
import { fitCurve } from './curves.js';
import { perp, cross, add, mul } from './fit.js';

/**
 * Devuelve { vertices, edges, faces, shells }.
 *  - vertices[i] = [x, y, z]
 *  - edges[i] = { v0, v1, curve, dev }
 *  - faces[i] = { S, sign, loops: [[{ edge, forward }]] , shell }
 *  - shells[i] = { faces: [índices], volume }
 */
export function buildBrep(mesh, seg, comp) {
  const { V, T, adj, nt } = mesh, { reg, prims, tol } = seg;
  const P = i => [V[i * 3], V[i * 3 + 1], V[i * 3 + 2]];
  const next = h => (h % 3 === 2 ? h - 2 : h + 1);
  const regOf = h => reg[(h / 3) | 0];
  const boundary = h => regOf(h) !== regOf(adj[h]);
  for (let h = 0; h < nt * 3; h++) if (adj[h] < 0) throw new Error('La malla no es cerrada.');

  // Vértices esquina: aquellos con un número de aristas frontera distinto de 2.
  const count = new Int32Array(mesh.nv);
  for (let h = 0; h < nt * 3; h++) if (boundary(h)) count[T[h]]++;

  const vertices = [], vertexId = new Map();
  const vertexOf = v => { let id = vertexId.get(v); if (id === undefined) { id = vertices.length; vertexId.set(v, id); vertices.push(P(v)); } return id; };
  const edges = [], chainOf = new Map(); // primera semiarista -> { edge, forward }
  const faces = prims.map((p, i) => ({ S: p.S, sign: p.sign, loops: [], shell: -1, region: i }));
  const visited = new Uint8Array(nt * 3);

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
    const closedChain = start < 0;
    if (closedChain) { start = 0; for (let i = 1; i < loop.length; i++) if (T[loop[i]] < T[loop[start]]) start = i; }
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
        const pts = vs.map(P), closed = vs[0] === vs[vs.length - 1];
        const { curve, dev } = fitCurve(pts, closed, prims[A].S, prims[B].S, tol);
        const id = edges.length;
        edges.push({ v0: vertexOf(vs[0]), v1: vertexOf(vs[vs.length - 1]), curve, dev, n: pts.length });
        ref = { edge: id, forward: true };
        chainOf.set(adj[chain[chain.length - 1]], { edge: id, forward: false });
      }
      refs.push(ref); i = j;
    }
    faces[A].loops.push(refs);
  }

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
    } else throw new Error('Superficie sin contorno no admitida.');
  }
  void cross;
  return { vertices, edges, faces: faces.concat(extra), shells };
}
