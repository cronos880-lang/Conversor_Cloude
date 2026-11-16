// node tools/convert.mjs entrada.stl salida.step
import fs from 'fs'; import path from 'path';
import { parseSTL } from '../static/core/mesh.js';
import { convert } from '../static/core/convert.js';
const [, , input, output] = process.argv;
const buf = fs.readFileSync(input);
const t0 = performance.now();
const r = convert(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)), { name: path.basename(input).replace(/\.stl$/i, '') });
fs.writeFileSync(output, r.step);
console.log(JSON.stringify({ ...r.report, seconds: +((performance.now() - t0) / 1000).toFixed(2), bytes: r.step.length }));
