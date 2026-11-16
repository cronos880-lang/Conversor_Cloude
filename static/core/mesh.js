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
    {
      // Triángulos sin área (vértices alineados): se eliminan; la reparación cierra el hueco en T.
      const ux = coords[b * 3] - coords[a * 3], uy = coords[b * 3 + 1] - coords[a * 3 + 1], uz = coords[b * 3 + 2] - coords[a * 3 + 2];
      const vx = coords[c * 3] - coords[a * 3], vy = coords[c * 3 + 1] - coords[a * 3 + 1], vz = coords[c * 3 + 2] - coords[a * 3 + 2];
      const cr = Math.hypot(uy * vz - uz * vy, uz * vx - ux * vz, ux * vy - uy * vx);
      const l2 = Math.max(ux * ux + uy * uy + uz * uz, vx * vx + vy * vy + vz * vz);
      if (!(cr > 1e-9 * l2)) { dropped++; continue; }
    }
    const lo = Math.min(a, b, c), key = lo === a ? `${a},${b},${c}` : lo === b ? `${b},${c},${a}` : `${c},${a},${b}`;
    if (seen.has(key)) { dropped++; continue; }
    seen.add(key);
    tri[nt * 3] = a; tri[nt * 3 + 1] = b; tri[nt * 3 + 2] = c; nt++;
  }
  let mesh = finishMesh(new Float64Array(coords), tri.slice(0, nt * 3), dropped);
  if (mesh.open) mesh = repair(mesh);
  return mesh;
}

/** Bordes abiertos: lista de semiaristas sin gemela. */
function openHalfEdges(mesh) {
  const out = [];
  for (let h = 0; h < mesh.nt * 3; h++) if (mesh.adj[h] < 0 && !mesh.multi[h]) out.push(h);
  return out;
}

/**
 * Reparación ligera de mallas abiertas: suelda vértices casi coincidentes,
 * resuelve uniones en T y tapa huecos pequeños. Anota lo hecho en mesh.repairs.
 */
function repair(mesh) {
  const repairs = { welded: 0, split: 0, filled: 0 };
  let V = Array.from(mesh.V), T = Array.from(mesh.T);
  const rebuild = () => { const m = finishMesh(new Float64Array(V), Int32Array.from(T), mesh.dropped); V = Array.from(m.V); T = Array.from(m.T); return m; };
  const P = i => [V[i * 3], V[i * 3 + 1], V[i * 3 + 2]];
  const d = (a, b) => Math.hypot(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
  let m = mesh;
  // 1) Soldar vértices de borde casi coincidentes (tolerancias crecientes).
  for (const eps of [mesh.diag * 1e-6, mesh.diag * 1e-5, mesh.diag * 1e-4]) {
    const open = openHalfEdges(m); if (!open.length) break;
    const ids = [...new Set(open.flatMap(h => [m.T[h], m.T[((h / 3) | 0) * 3 + (h % 3 + 1) % 3]]))];
    const grid = new Map(), cell = eps * 2, key = (x, y, z) => `${x},${y},${z}`;
    const remap = new Map();
    for (const i of ids) {
      const p = P(i), c = p.map(x => Math.floor(x / cell)); let found = -1;
      for (let dx = -1; dx <= 1 && found < 0; dx++) for (let dy = -1; dy <= 1 && found < 0; dy++) for (let dz = -1; dz <= 1 && found < 0; dz++) {
        for (const j of grid.get(key(c[0] + dx, c[1] + dy, c[2] + dz)) || []) if (d(P(j), p) <= eps) { found = j; break; }
      }
      if (found >= 0) remap.set(i, found);
      else { const k = key(c[0], c[1], c[2]); (grid.get(k) || grid.set(k, []).get(k)).push(i); }
    }
    if (!remap.size) continue;
    T = T.map(v => (remap.has(v) ? remap.get(v) : v));
    repairs.welded += remap.size;
    // Quitar los triángulos que hayan quedado degenerados.
    const keep = []; for (let t = 0; t < T.length / 3; t++) { const a = T[t * 3], b = T[t * 3 + 1], c = T[t * 3 + 2]; if (a !== b && b !== c && a !== c) keep.push(a, b, c); }
    T = keep; m = rebuild();
  }
  // 2) Uniones en T: un vértice de borde situado sobre una arista de borde la divide.
  for (let pass = 0; pass < 20; pass++) {
    const open = openHalfEdges(m); if (!open.length) break;
    const verts = [...new Set(open.map(h => m.T[h]))];
    const eps = mesh.diag * 1e-5; let done = false;
    for (const h of open) {
      const t = (h / 3) | 0, k = h % 3, a = T[t * 3 + k], b = T[t * 3 + (k + 1) % 3], c = T[t * 3 + (k + 2) % 3];
      const pa = P(a), pb = P(b), len = d(pa, pb); let best = -1, bs = 0;
      for (const v of verts) {
        if (v === a || v === b) continue;
        const p = P(v), s = ((p[0] - pa[0]) * (pb[0] - pa[0]) + (p[1] - pa[1]) * (pb[1] - pa[1]) + (p[2] - pa[2]) * (pb[2] - pa[2])) / (len * len);
        if (s <= 1e-9 || s >= 1 - 1e-9) continue;
        const q = [pa[0] + s * (pb[0] - pa[0]), pa[1] + s * (pb[1] - pa[1]), pa[2] + s * (pb[2] - pa[2])];
        if (d(p, q) <= eps && (best < 0 || s < bs)) { best = v; bs = s; }
      }
      if (best < 0) continue;
      T[t * 3] = a; T[t * 3 + 1] = best; T[t * 3 + 2] = c; T.push(best, b, c);
      repairs.split++; done = true; break;
    }
    if (!done) break;
    m = rebuild();
  }
  // 3) Tapar huecos: cada bucle de borde se triangula recortando orejas.
  for (let pass = 0; pass < 4; pass++) {
    const open = openHalfEdges(m); if (!open.length) break;
    const nextOf = new Map(); let simple = true;
    for (const h of open) { const a = m.T[h], b = m.T[((h / 3) | 0) * 3 + (h % 3 + 1) % 3]; if (nextOf.has(b)) simple = false; nextOf.set(b, a); } // sentido inverso al de la malla
    if (!simple) break;
    const used = new Set(); let added = 0;
    for (const start of nextOf.keys()) {
      if (used.has(start)) continue;
      const loop = []; let v = start;
      while (v !== undefined && !used.has(v)) { used.add(v); loop.push(v); v = nextOf.get(v); }
      if (v !== start || loop.length < 3 || loop.length > 500) continue;
      // Normal del bucle (método de Newell) para decidir qué vértices son convexos.
      const n = [0, 0, 0];
      for (let i = 0; i < loop.length; i++) { const p = P(loop[i]), q = P(loop[(i + 1) % loop.length]); n[0] += (p[1] - q[1]) * (p[2] + q[2]); n[1] += (p[2] - q[2]) * (p[0] + q[0]); n[2] += (p[0] - q[0]) * (p[1] + q[1]); }
      const ring = loop.slice(); let guard = ring.length * ring.length + 10;
      while (ring.length > 3 && guard-- > 0) {
        let cut = -1, bestScore = -Infinity;
        for (let i = 0; i < ring.length; i++) {
          const a = P(ring[(i + ring.length - 1) % ring.length]), b = P(ring[i]), c = P(ring[(i + 1) % ring.length]);
          const u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]], w = [c[0] - b[0], c[1] - b[1], c[2] - b[2]];
          const cr = [u[1] * w[2] - u[2] * w[1], u[2] * w[0] - u[0] * w[2], u[0] * w[1] - u[1] * w[0]];
          const convex = cr[0] * n[0] + cr[1] * n[1] + cr[2] * n[2];
          if (convex <= 0) continue;
          // Oreja válida: ningún otro vértice del anillo dentro del triángulo.
          let inside = false;
          for (let j = 0; j < ring.length && !inside; j++) {
            if (j === i || j === (i + 1) % ring.length || j === (i + ring.length - 1) % ring.length) continue;
            const p = P(ring[j]); let all = true;
            for (const [x, y] of [[a, b], [b, c], [c, a]]) {
              const e = [y[0] - x[0], y[1] - x[1], y[2] - x[2]], f = [p[0] - x[0], p[1] - x[1], p[2] - x[2]];
              const z = [e[1] * f[2] - e[2] * f[1], e[2] * f[0] - e[0] * f[2], e[0] * f[1] - e[1] * f[0]];
              if (z[0] * n[0] + z[1] * n[1] + z[2] * n[2] < 0) { all = false; break; }
            }
            inside = all;
          }
          if (inside) continue;
          const score = convex / (d(a, b) * d(b, c) * d(c, a) + 1e-300);
          if (score > bestScore) { bestScore = score; cut = i; }
        }
        if (cut < 0) break;
        const L = ring.length;
        T.push(ring[(cut + L - 1) % L], ring[cut], ring[(cut + 1) % L]); added++;
        ring.splice(cut, 1);
      }
      if (ring.length === 3) { T.push(ring[0], ring[1], ring[2]); added++; }
    }
    if (!added) break;
    repairs.filled += added; m = rebuild();
  }
  m.repairs = repairs;
  return m;
}

export function finishMesh(V, T, dropped = 0) {
  const nv = V.length / 3, nt = T.length / 3;
  // Gemelos de semiarista: adj[3t+k] = semiarista opuesta a la arista (k -> k+1) o -1.
  const adj = new Int32Array(nt * 3).fill(-1);
  const edgeMap = new Map();
  let open = 0, nonManifold = 0, flipped = 0, shared = 0;
  const multi = new Uint8Array(nt * 3);
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
    } else if (list.length % 2 === 0) {
      // Arista compartida por varios sólidos: emparejar caras consecutivas alrededor de la arista.
      const a = T[list[0]], b = T[((list[0] / 3) | 0) * 3 + (list[0] % 3 + 1) % 3];
      const e = [V[b * 3] - V[a * 3], V[b * 3 + 1] - V[a * 3 + 1], V[b * 3 + 2] - V[a * 3 + 2]];
      const el = Math.hypot(e[0], e[1], e[2]) || 1; e[0] /= el; e[1] /= el; e[2] /= el;
      const ref = Math.abs(e[0]) < 0.6 ? [1, 0, 0] : Math.abs(e[1]) < 0.6 ? [0, 1, 0] : [0, 0, 1];
      const x = [e[1] * ref[2] - e[2] * ref[1], e[2] * ref[0] - e[0] * ref[2], e[0] * ref[1] - e[1] * ref[0]];
      const xl = Math.hypot(x[0], x[1], x[2]); x[0] /= xl; x[1] /= xl; x[2] /= xl;
      const y = [e[1] * x[2] - e[2] * x[1], e[2] * x[0] - e[0] * x[2], e[0] * x[1] - e[1] * x[0]];
      const items = list.map(h => {
        const t = (h / 3) | 0, c = T[t * 3 + (h % 3 + 2) % 3];
        const r = [V[c * 3] - V[a * 3], V[c * 3 + 1] - V[a * 3 + 1], V[c * 3 + 2] - V[a * 3 + 2]];
        return { h, plus: T[h] === a, angle: Math.atan2(r[0] * y[0] + r[1] * y[1] + r[2] * y[2], r[0] * x[0] + r[1] * x[1] + r[2] * x[2]) };
      }).sort((p, q) => p.angle - q.angle);
      // Una cara "+" tiene el material hacia ángulos decrecientes: su pareja es la "-" anterior.
      let okPairs = true; const n = items.length;
      for (let i = 0; i < n && okPairs; i++) if (items[i].plus) {
        const prev = items[(i + n - 1) % n];
        if (prev.plus) okPairs = false; else { adj[items[i].h] = prev.h; adj[prev.h] = items[i].h; }
      }
      if (!okPairs) {
        // Caras coincidentes o casi: emparejar cada "+" con la "-" libre más cercana en ángulo.
        const plus = items.filter(it => it.plus), minus = items.filter(it => !it.plus);
        if (plus.length === minus.length) {
          const taken = new Set(); okPairs = true;
          for (const p of plus) {
            let best = null, gap = Infinity;
            for (const q of minus) { if (taken.has(q)) continue; let g = p.angle - q.angle; if (g < 0) g += 2 * Math.PI; if (g < gap) { gap = g; best = q; } }
            taken.add(best); adj[p.h] = best.h; adj[best.h] = p.h;
          }
        }
      }
      if (!okPairs) { for (const it of items) { adj[it.h] = -1; multi[it.h] = 1; } nonManifold++; } else shared++;
    } else { for (const h of list) multi[h] = 1; nonManifold++; }
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
  return { V, T, N, A, adj, multi, nv, nt, open, nonManifold, flipped, shared, zeroArea, dropped, volume, area, min, max, diag };
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
