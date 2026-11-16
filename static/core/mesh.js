// Malla indexada: soldadura de vértices, adyacencias y comprobaciones básicas.

/** Lee un STL binario o ASCII y devuelve Float32Array con 9 números por triángulo. */
export function parseSTL(buffer) {
  const bytes = new Uint8Array(buffer);
  if (bytes.length >= 84) {
    const count = new DataView(buffer).getUint32(80, true);
    if (84 + count * 50 === bytes.length) {
      const out = new Float32Array(count * 9), view = new DataView(buffer);
      for (let i = 0; i < count; i++) {
        const base = 84 + i * 50 + 12;
        for (let k = 0; k < 9; k++) out[i * 9 + k] = view.getFloat32(base + k * 4, true);
      }
      return out;
    }
  }
  const text = new TextDecoder().decode(bytes), out = [];
  const re = /vertex\s+(\S+)\s+(\S+)\s+(\S+)/g;
  let m;
  while ((m = re.exec(text))) out.push(+m[1], +m[2], +m[3]);
  return new Float32Array(out);
}

/**
 * Construye la malla indexada a partir de posiciones (9 números por triángulo).
 * Suelda vértices idénticos y elimina triángulos degenerados o duplicados.
 */
export function buildMesh(positions, scale = 1) {
  const nt0 = positions.length / 9;
  const map = new Map(), coords = [];
  const tri = new Int32Array(nt0 * 3);
  const f32 = new Float32Array(3), u32 = new Uint32Array(f32.buffer);
  for (let i = 0; i < nt0 * 3; i++) {
    f32[0] = positions[i * 3] + 0; f32[1] = positions[i * 3 + 1] + 0; f32[2] = positions[i * 3 + 2] + 0;
    if (f32[0] === 0) f32[0] = 0; if (f32[1] === 0) f32[1] = 0; if (f32[2] === 0) f32[2] = 0;
    const key = u32[0] + ',' + u32[1] + ',' + u32[2];
    let id = map.get(key);
    if (id === undefined) {
      id = coords.length / 3; map.set(key, id);
      coords.push(f32[0] * scale, f32[1] * scale, f32[2] * scale);
    }
    tri[i] = id;
  }
  // Quitar degenerados y duplicados exactos.
  const seen = new Set(); let nt = 0, dropped = 0;
  for (let t = 0; t < nt0; t++) {
    const a = tri[t * 3], b = tri[t * 3 + 1], c = tri[t * 3 + 2];
    if (a === b || b === c || a === c) { dropped++; continue; }
    const lo = Math.min(a, b, c), key = lo === a ? `${a},${b},${c}` : lo === b ? `${b},${c},${a}` : `${c},${a},${b}`;
    if (seen.has(key)) { dropped++; continue; }
    seen.add(key);
    tri[nt * 3] = a; tri[nt * 3 + 1] = b; tri[nt * 3 + 2] = c; nt++;
  }
  return finishMesh(new Float64Array(coords), tri.slice(0, nt * 3), dropped);
}

export function finishMesh(V, T, dropped = 0) {
  const nv = V.length / 3, nt = T.length / 3;
  // Gemelos de semiarista: adj[3t+k] = semiarista opuesta a la arista (k -> k+1) o -1.
  const adj = new Int32Array(nt * 3).fill(-1);
  const edgeMap = new Map();
  let open = 0, nonManifold = 0, flipped = 0;
  for (let h = 0; h < nt * 3; h++) {
    const t = (h / 3) | 0, k = h % 3;
    const a = T[h], b = T[t * 3 + (k + 1) % 3];
    const key = a < b ? a * nv + b : b * nv + a;
    const list = edgeMap.get(key);
    if (list === undefined) edgeMap.set(key, [h]); else list.push(h);
  }
  for (const list of edgeMap.values()) {
    if (list.length === 1) open++;
    else if (list.length === 2) {
      const [h1, h2] = list;
      if (T[h1] === T[h2]) flipped++;
      adj[h1] = h2; adj[h2] = h1;
    } else nonManifold++;
  }
  // Normales y áreas.
  const N = new Float64Array(nt * 3), A = new Float64Array(nt);
  let volume = 0, area = 0;
  const min = [Infinity, Infinity, Infinity], max = [-Infinity, -Infinity, -Infinity];
  for (let i = 0; i < nv; i++) for (let j = 0; j < 3; j++) {
    const v = V[i * 3 + j]; if (v < min[j]) min[j] = v; if (v > max[j]) max[j] = v;
  }
  const o = [(min[0] + max[0]) / 2, (min[1] + max[1]) / 2, (min[2] + max[2]) / 2];
  let zeroArea = 0;
  for (let t = 0; t < nt; t++) {
    const a = T[t * 3] * 3, b = T[t * 3 + 1] * 3, c = T[t * 3 + 2] * 3;
    const ux = V[b] - V[a], uy = V[b + 1] - V[a + 1], uz = V[b + 2] - V[a + 2];
    const vx = V[c] - V[a], vy = V[c + 1] - V[a + 1], vz = V[c + 2] - V[a + 2];
    const nx = uy * vz - uz * vy, ny = uz * vx - ux * vz, nz = ux * vy - uy * vx;
    const len = Math.hypot(nx, ny, nz);
    A[t] = len / 2; area += len / 2;
    if (len > 0) { N[t * 3] = nx / len; N[t * 3 + 1] = ny / len; N[t * 3 + 2] = nz / len; } else zeroArea++;
    volume += ((V[a] - o[0]) * nx + (V[a + 1] - o[1]) * ny + (V[a + 2] - o[2]) * nz) / 6;
  }
  const diag = Math.hypot(max[0] - min[0], max[1] - min[1], max[2] - min[2]);
  return { V, T, N, A, adj, nv, nt, open, nonManifold, flipped, zeroArea, dropped, volume, area, min, max, diag };
}

/** Etiqueta de componente conexa por triángulo. */
export function components(mesh) {
  const { nt, adj } = mesh, comp = new Int32Array(nt).fill(-1);
  let count = 0; const stack = [];
  for (let s = 0; s < nt; s++) {
    if (comp[s] >= 0) continue;
    comp[s] = count; stack.push(s);
    while (stack.length) {
      const t = stack.pop();
      for (let k = 0; k < 3; k++) {
        const h = adj[t * 3 + k]; if (h < 0) continue;
        const u = (h / 3) | 0; if (comp[u] < 0) { comp[u] = count; stack.push(u); }
      }
    }
    count++;
  }
  return { comp, count };
}
