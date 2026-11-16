// Prueba del flujo real en navegador: node tests/browser_real.mjs carpeta_stl carpeta_salida [filtro]
// Requiere un servidor estático en http://127.0.0.1:8011 sirviendo la raíz del proyecto.
import fs from 'fs'; import path from 'path'; import { createRequire } from 'module';
const require = createRequire(import.meta.url); const { chromium } = require('playwright');
const [, , dir, out, only] = process.argv; fs.mkdirSync(out, { recursive: true });
const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
const page = await browser.newPage(); const errors = [];
page.on('pageerror', e => errors.push(String(e)));
await page.goto(process.env.URL || 'http://127.0.0.1:8011/');
let ok = 0, total = 0;
for (const f of fs.readdirSync(dir).filter(f => /\.stl$/i.test(f) && (!only || f.includes(only))).sort()) {
  total++;
  await page.locator('#file').setInputFiles(path.join(dir, f));
  await page.waitForFunction("!document.querySelector('#convert').disabled || /no válido|Máximo|dimensiones/.test(document.querySelector('#status').textContent)", null, { timeout: 120000 });
  const t0 = Date.now();
  const download = page.waitForEvent('download', { timeout: 900000 }).catch(() => null);
  await page.locator('#convert').click();
  const failed = page.waitForFunction("document.querySelector('#cancel').hidden && !document.querySelector('#status').textContent.includes('STEP generado')", null, { timeout: 900000 }).then(() => null).catch(() => null);
  const done = await Promise.race([download, failed]);
  const status = (await page.locator('#status').innerText()).replace(/\n/g, ' | ');
  if (done) { await done.saveAs(path.join(out, f.replace(/\.stl$/i, '.step'))); ok++; }
  console.log(`${done ? 'OK ' : 'MAL'} ${f.slice(0, 34).padEnd(34)} ${((Date.now() - t0) / 1000).toFixed(1)}s  ${status}`);
}
console.log(`${ok}/${total} descargas; errores de página: ${errors.length ? errors.join('; ') : 'ninguno'}`);
await browser.close();
