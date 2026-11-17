// Dibuja la segmentación de un STL en un PNG: node tools/render.mjs pieza.stl salida.png [zoom] [tx,ty,tz]
import fs from 'fs'; import { createRequire } from 'module';
import { parseSTL, buildMesh } from '../static/core/mesh.js';
import { segment, defaultTolerance, mergeSameSurfaces } from '../static/core/segment.js';
import { findExtrusions } from '../static/core/extrude.js';
const require = createRequire(import.meta.url); const { chromium } = require('playwright');
const [, , file, out, zoom, target] = process.argv;
const buf = fs.readFileSync(file);
const m = buildMesh(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
const sg = segment(m, process.env.TOLK ? { tol: defaultTolerance(m) * +process.env.TOLK } : {}); mergeSameSurfaces(m, sg); findExtrusions(m, sg); const { reg, prims } = sg;
const hue = { extrusion: 0.15, plane: null, cylinder: 0.6, cone: 0.08, sphere: 0.33, torus: 0.85 };
const hsl = (h, s, l) => { const a = s * Math.min(l, 1 - l), f = n => { const k = (n + h * 12) % 12; return l - a * Math.max(-1, Math.min(k - 3, 9 - k, 1)); }; return [f(0), f(8), f(4)]; };
const colors = prims.map((p, i) => { const j = ((i * 2654435761) >>> 0) / 4294967296; return hue[p.S.type] === null ? hsl(0, 0, 0.6 + 0.3 * j) : hsl(hue[p.S.type] + (j - 0.5) * 0.1, 0.9, 0.4 + 0.3 * j); });
const pos = [], col = [], lines = [];
for (let t = 0; t < m.nt; t++) for (let k = 0; k < 3; k++) {
  const v = m.T[t * 3 + k] * 3; pos.push(m.V[v], m.V[v + 1], m.V[v + 2]); col.push(...colors[reg[t]]);
  const h = m.adj[t * 3 + k];
  if (h < 0 || (reg[(h / 3) | 0] !== reg[t] && h > t * 3 + k)) { const w = m.T[t * 3 + (k + 1) % 3] * 3; lines.push(m.V[v], m.V[v + 1], m.V[v + 2], m.V[w], m.V[w + 1], m.V[w + 2]); }
}
const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader', '--use-gl=angle', '--use-angle=swiftshader'] });
const page = await browser.newPage({ viewport: { width: 1500, height: 720 } });
await page.goto('http://127.0.0.1:8011/tools/render.html'); await page.waitForFunction('window.draw');
await page.evaluate(d => window.draw(d, 750, 720), { pos, col, lines, zoom: +zoom || 1, target: target ? target.split(',').map(Number) : null });
await page.screenshot({ path: out }); await browser.close();
console.log(out, prims.length, 'regiones');
