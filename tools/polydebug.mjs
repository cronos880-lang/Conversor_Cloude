import fs from 'fs';
import { parseSTL } from '../static/core/mesh.js';
import { convert } from '../static/core/convert.js';
const buf = fs.readFileSync(process.argv[2]);
const r = convert(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
const stats = {}; let shown = 0;
for (const c of r.brep.chainInfo) { const k = `${c.A}/${c.B} ${c.kind}`; stats[k] = (stats[k] || 0) + 1; if (c.kind === 'poly' && shown++ < 8) console.log(JSON.stringify(c)); }
console.log(stats);
