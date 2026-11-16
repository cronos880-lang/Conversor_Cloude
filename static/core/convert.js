// Conversión completa: posiciones STL -> texto STEP con superficies analíticas.
import { buildMesh, components } from './mesh.js';
import { segment, demote } from './segment.js';
import { buildBrep } from './brep.js';
import { writeStep } from './step.js';
import { validate } from './validate.js';

function prepare(positions, scale, onStatus) {
  onStatus('Analizando la malla…');
  const mesh = buildMesh(positions, scale);
  if (!mesh.nt) throw new Error('STL vacío o no válido.');
  if (mesh.open || mesh.nonManifold) throw new Error('Se requiere una malla cerrada. Repara el STL antes de convertir.');
  if (mesh.flipped) throw new Error('La malla tiene triángulos con orientación incoherente.');
  if (!(mesh.volume > 1e-9 * mesh.diag ** 3)) throw new Error('El STL no encierra un volumen válido.');
  const comp = components(mesh);
  onStatus('Reconociendo superficies…');
  const seg = segment(mesh);
  return { mesh, comp, seg };
}

function build(state, name, polyFaces) {
  const { mesh, comp, seg } = state;
  const brep = buildBrep(mesh, seg, comp, { polyFaces });
  let dev = 0; for (const e of brep.edges) if (!e.poly) dev = Math.max(dev, e.dev);
  const step = writeStep(brep, name, Math.max(1e-6, Math.min(dev, 1e-2)));
  const types = {}; for (const f of brep.faces) if (f.loops.length) types[f.S.type] = (types[f.S.type] || 0) + 1;
  const curves = {}; for (const e of brep.edges) curves[e.poly ? 'poly' : e.curve.type] = (curves[e.poly ? 'poly' : e.curve.type] || 0) + 1;
  const faces = brep.faces.filter(f => f.loops.length).length;
  return { step, brep, report: { triangles: mesh.nt, faces, types, curves, edgeDeviation: dev, volume: mesh.volume, area: mesh.area, tolerance: seg.tol } };
}

/** Conversión sin validar (para pruebas y diagnóstico). */
export function convert(positions, { scale = 1, name = 'pieza', onStatus = () => {} } = {}) {
  const state = prepare(positions, scale, onStatus);
  onStatus('Construyendo el sólido…');
  return { ...build(state, name, new Set()), ...state };
}

/**
 * Conversión validada con Open CASCADE. Si alguna cara no es válida se degrada
 * primero a contornos poligonales y después a facetas planas, hasta obtener un
 * sólido válido que conserva el volumen.
 */
export function convertValidated(positions, oc, { scale = 1, name = 'pieza', onStatus = () => {}, trace = null } = {}) {
  const state = prepare(positions, scale, onStatus);
  const { mesh, seg } = state;
  const polyFaces = new Set(); let demoted = 0, rounds = 0, good = null;
  for (; rounds < 7; rounds++) {
    onStatus(rounds ? `Corrigiendo caras problemáticas (intento ${rounds})…` : 'Construyendo el sólido…');
    const t0 = performance.now();
    const out = build(state, name, polyFaces);
    const regionArea = new Float64Array(seg.prims.length);
    for (let t = 0; t < mesh.nt; t++) regionArea[seg.reg[t]] += mesh.A[t];
    onStatus('Comprobando el sólido…');
    const check = validate(oc, out.step, out.brep, regionArea);
    const volumeOk = Math.abs(check.volume - mesh.volume) <= 0.02 * Math.abs(mesh.volume);
    const areaOk = Math.abs(check.area - mesh.area) <= 0.02 * mesh.area;
    const solid = check.valid && volumeOk && areaOk;
    if (trace) trace({ round: rounds, ms: Math.round(performance.now() - t0), valid: check.valid, readable: check.readable, faces: check.faces, mine: out.report.faces, unmatched: check.unmatched, free: check.free, bad: [...check.bad].slice(0, 12), nbad: check.bad.size, vol: +(check.volume / mesh.volume).toFixed(4), area: +(check.area / mesh.area).toFixed(4), error: check.error, why: (check.why || []).slice(0, 6).map(w => w.concat(seg.prims[w[0]].S.type)) });
    if (solid && (!good || check.bad.size <= good.check.bad.size)) good = { ...out, check, rounds };
    if (solid && !check.bad.size) break;
    // Caras a degradar: las señaladas o, si el fallo es global, todas las que aún no son facetas.
    let bad = [...check.bad];
    if (!bad.length) bad = seg.prims.map((p, i) => i).filter(i => !(polyFaces.has(i) && seg.prims[i].facet));
    const toDemote = new Set(); let progress = false, neighbours = null;
    for (const r of bad) {
      const gross = check.gross && check.gross.has(r) && !seg.prims[r].facet && seg.prims[r].S.type !== 'plane';
      if (gross) { polyFaces.add(r); toDemote.add(r); }                                    // área disparatada: facetas ya
      else if (!polyFaces.has(r)) { polyFaces.add(r); progress = true; }                  // 1.º contorno poligonal
      else if (!seg.prims[r].facet) toDemote.add(r);   // 2.º facetas planas
      else {
        // Una faceta que sigue fallando: degradar las caras vecinas que aún no son facetas.
        if (!neighbours) {
          neighbours = new Map();
          for (let h = 0; h < mesh.nt * 3; h++) {
            const a = seg.reg[(h / 3) | 0], b = seg.reg[(mesh.adj[h] / 3) | 0];
            if (a !== b) (neighbours.get(a) || neighbours.set(a, new Set()).get(a)).add(b);
          }
        }
        for (const q of neighbours.get(r) || []) if (!seg.prims[q].facet) { if (!polyFaces.has(q)) polyFaces.add(q); else toDemote.add(q); progress = true; }
      }
    }
    if (toDemote.size) {
      progress = true; demoted += demote(mesh, seg, toDemote);
      for (let i = 0; i < seg.prims.length; i++) if (seg.prims[i].facet) polyFaces.add(i);
    }
    if (!progress) break;
  }
  if (!good) {
    // Último recurso: todo el sólido con facetas planas (como una conversión facetada clásica).
    onStatus('Usando caras planas en toda la pieza…');
    const all = new Set(seg.prims.map((p, i) => i));
    demoted += demote(mesh, seg, all);
    const every = new Set(seg.prims.map((p, i) => i));
    const out = build(state, name, every);
    const regionArea = new Float64Array(seg.prims.length);
    for (let t = 0; t < mesh.nt; t++) regionArea[seg.reg[t]] += mesh.A[t];
    const check = validate(oc, out.step, out.brep, regionArea);
    if (trace) trace({ round: 'facetas', valid: check.valid, vol: +(check.volume / mesh.volume).toFixed(4), area: +(check.area / mesh.area).toFixed(4), nbad: check.bad.size, error: check.error });
    if (check.valid && Math.abs(check.volume - mesh.volume) <= 0.02 * Math.abs(mesh.volume)) good = { ...out, check, rounds };
  }
  if (!good) throw new Error('No se pudo obtener un sólido válido a partir de este STL.');
  const last = good; rounds = good.rounds;
  last.report.rounds = rounds; last.report.demoted = demoted;
  if (mesh.repairs) last.report.repairs = mesh.repairs;
  last.report.checkedVolume = last.check.volume;
  return { ...last, ...state };
}
