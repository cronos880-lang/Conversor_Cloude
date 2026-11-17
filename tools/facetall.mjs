// STEP completamente facetado (sin reconocer superficies): node tools/facetall.mjs pieza.stl salida.step
import fs from 'fs';
import { parseSTL, buildMesh, components } from '../static/core/mesh.js';
import { segment, demote } from '../static/core/segment.js';
import { buildBrep } from '../static/core/brep.js';
import { writeStep } from '../static/core/step.js';
const buf = fs.readFileSync(process.argv[2]);
const mesh = buildMesh(parseSTL(buf.buffer.slice(buf.byteOffset, buf.byteOffset + buf.length)));
const seg = process.env.FROMSEG ? segment(mesh) : { reg: new Int32Array(mesh.nt), prims: [{ S: { type: 'plane', n: [0, 0, 1], d: 0 }, sign: 1 }], tol: 0.01 };
demote(mesh, seg, new Set(seg.prims.map((p, i) => i)));
const brep = buildBrep(mesh, seg, components(mesh), { polyFaces: new Set(seg.prims.map((p, i) => i)) });
fs.writeFileSync(process.argv[3], writeStep(brep, 'facetas', 1e-6));
console.log('triángulos', mesh.nt, 'caras', brep.faces.filter(f => f.loops.length).length, JSON.stringify(mesh.repairs || {}));
