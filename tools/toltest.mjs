import fs from 'fs'; import path from 'path';
import { parseSTL, buildMesh } from '../static/core/mesh.js';
import { segment, defaultTolerance } from '../static/core/segment.js';
const dir = process.argv[2], names = process.argv.slice(3);
for (const f of fs.readdirSync(dir).filter(f => /\.stl$/i.test(f) && names.some(n => f.includes(n))).sort()) {
  const buf = fs.readFileSync(path.join(dir, f));
  const m = buildMesh(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
  const base = defaultTolerance(m); const out = [];
  for (const k of (process.env.KS || "1,3,10,30,100").split(",").map(Number)) { const t0 = performance.now(); const { prims } = segment(m, { tol: base * k }); const c = {}; for (const p of prims) c[p.S.type[0] + p.S.type[1]] = (c[p.S.type[0] + p.S.type[1]] || 0) + 1; out.push(`x${k}: ${prims.length} ${JSON.stringify(c).replace(/"/g, '')} ${((performance.now() - t0) / 1000).toFixed(1)}s`); }
  console.log(f.slice(0, 30), 'base', base.toExponential(1), '\n   ' + out.join('\n   '));
}
