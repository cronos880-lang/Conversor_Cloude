// Conversión completa: posiciones STL -> texto STEP con superficies analíticas.
import { buildMesh, components } from './mesh.js';
import { segment } from './segment.js';
import { buildBrep } from './brep.js';
import { writeStep } from './step.js';

export function convert(positions, { scale = 1, name = 'pieza', onStatus = () => {} } = {}) {
  onStatus('Analizando la malla…');
  const mesh = buildMesh(positions, scale);
  if (!mesh.nt) throw new Error('STL vacío o no válido.');
  if (mesh.open || mesh.nonManifold) throw new Error('Se requiere una malla cerrada. Repara el STL antes de convertir.');
  if (mesh.flipped) throw new Error('La malla tiene triángulos con orientación incoherente.');
  const comp = components(mesh);
  onStatus('Reconociendo superficies…');
  const seg = segment(mesh);
  onStatus('Construyendo el sólido…');
  const brep = buildBrep(mesh, seg, comp);
  let dev = 0; for (const e of brep.edges) dev = Math.max(dev, e.dev);
  const step = writeStep(brep, name, Math.max(1e-6, Math.min(seg.tol, 1e-3)));
  const types = {}; for (const f of brep.faces) types[f.S.type] = (types[f.S.type] || 0) + 1;
  const curves = {}; for (const e of brep.edges) curves[e.curve.type] = (curves[e.curve.type] || 0) + 1;
  return { step, mesh, seg, brep, report: { triangles: mesh.nt, faces: brep.faces.length, types, curves, edgeDeviation: dev, volume: mesh.volume, area: mesh.area, tolerance: seg.tol } };
}
