// Piezas sintéticas por revolución para pruebas controladas.
// node tools/synth.mjs salida.stl segmentos [ruido_mm]
import fs from 'fs';
export function revolve(profile, n, noise = 0) {
  // profile: [[r, z], ...] de un eje a otro (r = 0 en los extremos).
  const tris = [], rnd = () => (Math.random() - 0.5) * 2 * noise;
  const P = (r, z, i) => [r * Math.cos(2 * Math.PI * i / n) + rnd(), r * Math.sin(2 * Math.PI * i / n) + rnd(), z + rnd()];
  const pts = profile.map(([r, z]) => Array.from({ length: n }, (_, i) => r === 0 ? [0, 0, z] : P(r, z, i)));
  for (let k = 0; k + 1 < profile.length; k++) for (let i = 0; i < n; i++) {
    const j = (i + 1) % n, a = pts[k][i], b = pts[k][j], c = pts[k + 1][j], d = pts[k + 1][i];
    if (profile[k][0] !== 0) tris.push([a, c, b]);
    if (profile[k + 1][0] !== 0) tris.push([a, d, c]);
  }
  return tris;
}
export function toSTL(tris) {
  const b = Buffer.alloc(84 + tris.length * 50); b.writeUInt32LE(tris.length, 80);
  tris.forEach((t, i) => { let o = 84 + i * 50 + 12; for (const p of t) for (const v of p) { b.writeFloatLE(v, o); o += 4; } });
  return b;
}
// Placa con tetón y redondeo en la base: plano, cilindro exterior de la placa, plano superior, toro, cilindro del tetón, plano.
export function boss(fillet = 1.5, arcSeg = 6) {
  const prof = [[0, 0], [20, 0], [20, 2], [5 + fillet, 2]];
  for (let k = 1; k < arcSeg; k++) { const t = (Math.PI / 2) * k / arcSeg; prof.push([5 + fillet - fillet * Math.sin(t), 2 + fillet - fillet * Math.cos(t)]); }
  prof.push([5, 2 + fillet], [5, 8], [0, 8]);
  return prof;
}
if (process.argv[1].endsWith('synth.mjs')) {
  const [, , out, n, noise, arc] = process.argv;
  fs.writeFileSync(out, toSTL(revolve(boss(1.5, +(arc || 6)), +n, +(noise || 0))));
}
