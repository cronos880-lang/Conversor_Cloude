// Prueba del flujo real en navegador con la interfaz R21:
// node tests/browser_real.mjs carpeta_stl carpeta_salida [filtro]
// Requiere un servidor estático (por defecto http://127.0.0.1:8011/) sirviendo la raíz del proyecto.
import fs from 'fs'; import path from 'path'; import { createRequire } from 'module';
const require = createRequire(import.meta.url); const { chromium } = require('playwright');
const [, , dir, out, only] = process.argv; fs.mkdirSync(out, { recursive: true });
const browser = await chromium.launch({ args: ['--enable-unsafe-swiftshader'] });
const page = await browser.newPage({ viewport: { width: 1400, height: 800 } }); const errors = [];
page.on('pageerror', e => errors.push(String(e)));
page.on('dialog', d => d.accept(d.defaultValue()));
await page.goto(process.env.URL || 'http://127.0.0.1:8011/');
let ok = 0, total = 0;
for (const f of fs.readdirSync(dir).filter(f => /\.stl$/i.test(f) && (!only || f.includes(only))).sort()) {
  total++;
  await page.locator('#fileInput').setInputFiles(path.join(dir, f));
  await page.waitForFunction(() => !document.querySelector('#convertBtn').disabled, null, { timeout: 120000 });
  const t0 = Date.now();
  await page.locator('#convertBtn').click();
  await page.waitForFunction(() => !document.querySelector('#exportBtn').disabled || document.querySelector('#conversionProgress').dataset.state === 'failed', null, { timeout: 900000 });
  const seconds = ((Date.now() - t0) / 1000).toFixed(1);
  const details = await page.locator('#details').innerText();
  const done = !(await page.locator('#exportBtn').isDisabled());
  if (done) {
    const [download] = await Promise.all([page.waitForEvent('download'), page.locator('#exportBtn').click()]);
    await download.saveAs(path.join(out, f.replace(/\.stl$/i, '.step'))); ok++;
  }
  const toast = done ? '' : ' | ' + await page.locator('#error').innerText();
  console.log(`${done ? 'OK ' : 'MAL'} ${f.slice(0, 34).padEnd(34)} ${seconds}s  ${details}${toast}`);
}
if (process.env.SHOT) await page.screenshot({ path: process.env.SHOT });
console.log(`${ok}/${total} descargas; errores de página: ${errors.length ? errors.join('; ') : 'ninguno'}`);
await browser.close();
