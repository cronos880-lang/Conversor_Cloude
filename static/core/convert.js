// Conversión completa: posiciones STL -> texto STEP con superficies analíticas.
import { buildMesh, components, dropFlatComponents, finishMesh } from './mesh.js';
import { segment, demote, mergeSameSurfaces, splitRegion, defaultTolerance } from './segment.js';
import { buildBrep } from './brep.js';
import { writeStep } from './step.js';
import { validate, reexport } from './validate.js';
import { findExtrusions } from './extrude.js';
import { regularize } from './regularize.js';
import { exactFacetJoins } from './seam.js';

// Dividir caras que fallan: desactivado (en el chasis un trozo de toro se interpretaba al revés).
const SPLIT_FACES = !!globalThis.process?.env?.SPLIT_FACES;
const EDGE_LIMIT = Number(globalThis.process?.env?.EDGE_LIMIT) || 0.1; // mm (EDGE_LIMIT: solo pruebas)

/** Regiones curvas con alguna arista cuya desviación supera `limit`. */
function looseRegions(brep, prims, limit) {
  const out = new Set();
  for (const f of brep.faces) {
    if (f.S.type === 'plane' || prims[f.region].facet) continue;
    if (f.loops.some(loop => loop.some(ref => brep.edges[ref.edge].dev > limit))) out.add(f.region);
  }
  return out;
}

function prepare(positions, scale, onStatus, tolFactor = 1) {
  onStatus('Analizando la malla…');
  let mesh = buildMesh(positions, scale);
  if (!mesh.open && !mesh.nonManifold) mesh = dropFlatComponents(mesh);
  if (!mesh.nt) throw new Error('STL vacío o no válido.');
  if (mesh.open || mesh.nonManifold) throw new Error('Se requiere una malla cerrada. Repara el STL antes de convertir.');
  if (mesh.flipped) throw new Error('La malla tiene triángulos con orientación incoherente.');
  if (mesh.volume < 0) {
    // Malla con todas las normales hacia dentro: invertir todos los triángulos.
    const T = mesh.T.slice(); for (let t = 0; t < mesh.nt; t++) { const b = T[t * 3 + 1]; T[t * 3 + 1] = T[t * 3 + 2]; T[t * 3 + 2] = b; }
    const repairs = mesh.repairs; mesh = finishMesh(mesh.V, T, mesh.dropped); mesh.repairs = { ...(repairs || {}), flipped: true };
  }
  if (!(mesh.volume > 1e-9 * mesh.diag ** 3)) throw new Error('El STL no encierra un volumen válido.');
  const comp = components(mesh);
  onStatus('Reconociendo superficies…');
  const seg = segment(mesh, tolFactor !== 1 ? { tol: Math.min(0.06, defaultTolerance(mesh) * tolFactor) } : {});
  mergeSameSurfaces(mesh, seg);
  findExtrusions(mesh, seg);
  mergeSameSurfaces(mesh, seg);
  if (!globalThis.process?.env?.NOREG) regularize(mesh, seg);
  return { mesh, comp, seg };
}

function build(state, name, polyFaces) {
  const { mesh, comp, seg } = state;
  const moved = globalThis.process?.env?.NOJOIN ? null : exactFacetJoins(mesh, seg);
  const brep = buildBrep(mesh, seg, comp, { polyFaces, polyGap: state.polyGap, move: state.move, moved });
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
/**
 * `engine` es el motor Open CASCADE o una función que devuelve uno nuevo (promesa).
 * Con una función, en mallas grandes se usa un motor limpio en cada validación,
 * porque leer STEP enormes agota o corrompe la memoria del motor.
 */
export async function convertValidated(positions, engine, options = {}) {
  // Primer intento: vértices llevados con decisión a la intersección exacta de sus caras.
  const first = await attempt(positions, engine, options, 1e-9, 3);
  // Si hubo que degradar caras, probar con vértices más cerca de la malla y quitando también
  // los puntos casi repetidos de los bordes poligonales, y quedarse con el mejor.
  if (!first.report.demoted || first.mesh.nt > 50000) return first;
  try {
    const second = await attempt(positions, engine, { ...options, trace: null, dump: null }, 0.002, 1);
    return second.report.faces < first.report.faces ? second : first;
  } catch { return first; }
}

async function attempt(positions, engine, { scale = 1, name = 'pieza', onStatus = () => {}, trace = null, maxRounds = Infinity, dump = null, tolFactor = 1, facetRegions = null } = {}, polyGap = 1e-9, move = 1) {
  const factory = typeof engine === 'function' ? engine : null;
  let oc = factory ? await factory() : engine, fresh = true;
  const check_ = async (step, brep, regionArea) => {
    if (factory && !fresh) { onStatus('Preparando el motor CAD…'); oc = await factory(); }
    fresh = false;
    const r = validate(oc, step, brep, regionArea);
    // Tras un fallo interno el motor queda inutilizable: pedir uno nuevo la próxima vez.
    if (!factory || !(r.error || mesh_.nt > 50000)) fresh = true;
    return r;
  };
  let mesh_ = { nt: 0 };
  const state = prepare(positions, scale, onStatus, tolFactor); state.polyGap = polyGap; state.move = move;
  const { mesh, seg } = state; mesh_ = mesh;
  // Solo para diagnóstico: pasar a facetas desde el principio las regiones indicadas.
  if (facetRegions) { const sel = new Set(); seg.prims.forEach((p, i) => { if (facetRegions(p, i)) sel.add(i); }); if (sel.size) demote(mesh, seg, sel); }
  const polyFaces = new Set(); let demoted = 0, rounds = 0, good = null;
  // En mallas grandes cada ronda es cara: las caras defectuosas pasan directamente a facetas.
  const big = mesh.nt > 50000;
  for (; rounds < Math.min(maxRounds, big ? 10 : 7); rounds++) {
    onStatus(rounds ? `Corrigiendo caras problemáticas (intento ${rounds})…` : 'Construyendo el sólido…');
    const t0 = performance.now();
    const out = build(state, name, polyFaces);
    const regionArea = new Float64Array(seg.prims.length);
    for (let t = 0; t < mesh.nt; t++) regionArea[seg.reg[t]] += mesh.A[t];
    onStatus('Comprobando el sólido…');
    if (dump) dump(rounds, out.step);
    const check = await check_(out.step, out.brep, regionArea);
    // Holguras: una cara curva cuyos bordes se separan de la superficie más de lo que
    // admiten otros programas CAD (Onshape/Parasolid) se trata como defectuosa.
    for (const r of looseRegions(out.brep, seg.prims, EDGE_LIMIT)) { check.bad.add(r); if (!check.reason.has(r)) check.reason.set(r, 'holgura'); }
    const volumeOk = Math.abs(check.volume - mesh.volume) <= 0.02 * Math.abs(mesh.volume);
    const areaOk = Math.abs(check.area - mesh.area) <= 0.02 * mesh.area;
    const solid = check.valid && volumeOk && areaOk;
    if (trace) trace({ round: rounds, ms: Math.round(performance.now() - t0), valid: check.valid, readable: check.readable, faces: check.faces, mine: out.report.faces, unmatched: check.unmatched, free: check.free, bad: [...check.bad].slice(0, 12), nbad: check.bad.size, vol: +(check.volume / mesh.volume).toFixed(4), area: +(check.area / mesh.area).toFixed(4), error: check.error, why: (check.why || []).slice(0, 6).map(w => w.concat(seg.prims[w[0]].S.type)), kinds: (() => { const k = {}; for (const r of check.bad) { const key = (seg.prims[r].facet ? 'faceta' : seg.prims[r].S.type) + ':' + (check.reason.get(r) || '?'); k[key] = (k[key] || 0) + 1; } return k; })() });
    if (solid && (!good || check.bad.size <= good.check.bad.size)) good = { ...out, check, rounds };
    // Las facetas ya no se pueden degradar más: si solo quedan facetas señaladas en un sólido válido, se acepta.
    const pending = [...check.bad].filter(r => !seg.prims[r].facet);
    if (solid && !pending.length) { if (!good || good.check.bad.size) good = { ...out, check, rounds }; break; }
    // Caras a degradar: las señaladas o, si el fallo es global, todas las que aún no son facetas.
    let bad = [...check.bad];
    if (!bad.length) bad = seg.prims.map((p, i) => i).filter(i => !(polyFaces.has(i) && seg.prims[i].facet));
    const toDemote = new Set(); let progress = false, neighbours = null;
    const triCount = new Map(); for (let t = 0; t < mesh.nt; t++) triCount.set(seg.reg[t], (triCount.get(seg.reg[t]) || 0) + 1);
    const regionTris = r => triCount.get(r) || 0;
    for (const r of bad) {
      const gross = check.gross && check.gross.has(r) && !seg.prims[r].facet && seg.prims[r].S.type !== 'plane';
      // Redondeos pequeños (toros, esferas) defectuosos: facetas directamente. Las superficies
      // grandes solo pasan a facetas tras fallar también con contorno poligonal.
      const small = !seg.prims[r].facet && regionArea[r] < 0.01 * mesh.area && (seg.prims[r].S.type === 'torus' || seg.prims[r].S.type === 'sphere');
      const revol = ['torus', 'cylinder', 'cone'].includes(seg.prims[r].S.type) && !seg.prims[r].facet;
      if (revol && check.reason.get(r) === 'cara' && !seg.prims[r].seamTried) {
        // 0.º Costura de la superficie al lado opuesto de la cara (contorno que se cruzaba).
        seg.prims[r] = { ...seg.prims[r], seamAway: true, seamTried: true }; progress = true;
      }
      else if (SPLIT_FACES && revol && check.reason.get(r) === 'cara' && (seg.prims[r].split || 0) < 3 && regionTris(r) >= 40) {
        // 0.º bis: dividir la cara en partes más simples sobre la misma superficie.
        splitRegion(mesh, seg, r); progress = true;
      }
      else if (gross || small) { polyFaces.add(r); toDemote.add(r); }                                    // área disparatada: facetas ya
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
      if (trace) { const cnt = new Map(); for (let t = 0; t < mesh.nt; t++) if (toDemote.has(seg.reg[t])) cnt.set(seg.reg[t], (cnt.get(seg.reg[t]) || 0) + 1); trace({ degradadas: [...cnt].map(([r, n]) => `${seg.prims[r].S.type}:${n}tri:${(check.reason.get(r) || '')}`) }); }
      progress = true; demoted += demote(mesh, seg, toDemote);
      for (let i = 0; i < seg.prims.length; i++) if (seg.prims[i].facet) polyFaces.add(i);
    }
    if (!progress) break;
  }
  let fallbackGood = null;
  // Si casi toda la pieza quedó en facetas, las pocas caras curvas restantes dejan
  // holguras en sus bordes que Onshape marca como defectos: mejor todo en facetas exactas.
  if (good) {
    let facetArea = 0; for (let t = 0; t < mesh.nt; t++) if (seg.prims[seg.reg[t]].facet) facetArea += mesh.A[t];
    const curved = good.brep.faces.filter(f => f.loops.length && f.S.type !== 'plane').length;
    if (curved && facetArea > 0.8 * mesh.area) { fallbackGood = good; good = null; onStatus('Demasiadas zonas sin reconocer: usando caras planas exactas…'); }
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
    const check = await check_(out.step, out.brep, regionArea);
    if (trace) trace({ round: 'facetas', valid: check.valid, vol: +(check.volume / mesh.volume).toFixed(4), area: +(check.area / mesh.area).toFixed(4), nbad: check.bad.size, error: check.error });
    if (check.valid && Math.abs(check.volume - mesh.volume) <= 0.02 * Math.abs(mesh.volume)) good = { ...out, check, rounds };
  }
  // Si el facetado completo no sale válido, se conserva el resultado anterior, que sí lo era.
  if (!good) good = fallbackGood;
  if (!good) throw new Error('No se pudo obtener un sólido válido a partir de este STL.');
  const last = good; rounds = good.rounds;
  last.report.rounds = rounds; last.report.demoted = demoted;
  if (mesh.repairs) last.report.repairs = mesh.repairs;
  last.report.checkedVolume = last.check.volume;
  // STEP final escrito por el motor CAD (curvas de parámetros y costuras para otros CAD).
  onStatus('Escribiendo el STEP final…');
  const engineFinal = factory ? await factory() : oc;
  const final = globalThis.process?.env?.NOREEXPORT ? null : reexport(engineFinal, last.step);
  if (final) { last.step = final; last.report.reexported = true; }
  return { ...last, ...state };
}
