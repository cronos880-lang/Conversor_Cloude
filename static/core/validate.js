// Validación del STEP con Open CASCADE (WebAssembly): solidez, volumen y caras defectuosas.

const key = p => `${Math.round(p[0] * 1e5)},${Math.round(p[1] * 1e5)},${Math.round(p[2] * 1e5)}`;

/**
 * Lee el STEP con Open CASCADE y lo compara con el B-rep generado.
 * Devuelve { valid, volume, area, faces, bad } donde `bad` es el conjunto de
 * regiones cuya cara no es válida, falta o no tiene el área esperada.
 */
export function validate(oc, stepText, brep, regionArea) {
  const E = oc.TopAbs_ShapeEnum, done = oc.IFSelect_ReturnStatus.IFSelect_RetDone;
  oc.FS.writeFile('/check.step', stepText);
  const reader = new oc.STEPControl_Reader(), progress = new oc.Message_ProgressRange();
  const result = { valid: false, volume: 0, area: 0, faces: 0, bad: new Set(), readable: false };
  try {
    if (reader.ReadFile('/check.step') !== done || !reader.TransferRoots(progress)) return result;
    const shape = reader.OneShape();
    result.readable = true;
    const check = new oc.BRepCheck_Analyzer(shape, true, false, false);
    result.valid = check.IsValid(); check.delete();
    const props = new oc.GProp_GProps();
    oc.BRepGProp.VolumeProperties(shape, props, true, false, false); result.volume = props.Mass();
    oc.BRepGProp.SurfaceProperties(shape, props, false, false); result.area = props.Mass();
    props.delete();

    // Vértices del B-rep -> caras que los usan.
    const vertexFaces = new Map();
    brep.faces.forEach((f, i) => { for (const loop of f.loops) for (const ref of loop) for (const v of [brep.edges[ref.edge].v0, brep.edges[ref.edge].v1]) {
      const k = key(brep.vertices[v]); let set = vertexFaces.get(k); if (!set) vertexFaces.set(k, set = new Set()); set.add(i);
    } });
    const got = new Float64Array(brep.faces.length), seen = new Uint8Array(brep.faces.length);
    // Aristas libres: las que solo usa una cara indican que la envolvente no cierra.
    const edgeUse = new Map(); // hash -> [{ edge, count, faces }]
    const useEdge = (edge, faceIndex, reversedFace) => {
      if (oc.BRep_Tool.Degenerated(oc.TopoDS.Edge(edge))) return;
      const ends = []; const ev = new oc.TopExp_Explorer(edge, E.TopAbs_VERTEX, E.TopAbs_SHAPE);
      for (; ev.More(); ev.Next()) { const p = oc.BRep_Tool.Pnt(oc.TopoDS.Vertex(ev.Current())); ends.push(key([p.X(), p.Y(), p.Z()])); p.delete(); }
      ev.delete();
      const h = ends.sort().join('|'); let list = edgeUse.get(h); if (!list) edgeUse.set(h, list = []);
      const o = edge.Orientation().value; void reversedFace;
      for (const item of list) if (item.edge.IsSame(edge)) { item.count++; item.faces.add(faceIndex); item.orients.push(o); return; }
      list.push({ edge, count: 1, faces: new Set([faceIndex]), orients: [o] });
    };
    const faceRegion = [];
    const ex = new oc.TopExp_Explorer(shape, E.TopAbs_FACE, E.TopAbs_SHAPE);
    let index = 0; const unmatched = [];
    for (; ex.More(); ex.Next(), index++) {
      const face = oc.TopoDS.Face(ex.Current());
      const g = new oc.GProp_GProps(); oc.BRepGProp.SurfaceProperties(face, g, false, false);
      const area = g.Mass(); g.delete();
      const fc = new oc.BRepCheck_Analyzer(face, true, false, false), ok = fc.IsValid(); fc.delete();
      const ee = new oc.TopExp_Explorer(face, E.TopAbs_EDGE, E.TopAbs_SHAPE);
      const rev = face.Orientation().value === 1;
      for (; ee.More(); ee.Next()) useEdge(ee.Current(), index, rev);
      ee.delete();
      const votes = new Map(); let nv = 0;
      const vx = new oc.TopExp_Explorer(face, E.TopAbs_VERTEX, E.TopAbs_SHAPE);
      for (; vx.More(); vx.Next()) {
        const p = oc.BRep_Tool.Pnt(oc.TopoDS.Vertex(vx.Current())); nv++;
        const set = vertexFaces.get(key([p.X(), p.Y(), p.Z()])); p.delete();
        if (set) for (const i of set) votes.set(i, (votes.get(i) || 0) + 1);
      }
      vx.delete();
      let best = -1, score = 0, tie = false;
      for (const [i, c] of votes) { if (c > score) { best = i; score = c; tie = false; } else if (c === score) tie = true; }
      if (tie) {
        // Desempate: la cara candidata con área más parecida.
        let diff = Infinity;
        for (const [i, c] of votes) if (c === score) { const d = Math.abs(regionArea[brep.faces[i].region] - area) + (seen[i] ? 1e30 : 0); if (d < diff) { diff = d; best = i; } }
      }
      if (best < 0) { unmatched.push({ ok, area }); continue; }
      faceRegion[index] = brep.faces[best].region;
      seen[best] = 1; got[best] += area;
      if (!ok) result.bad.add(brep.faces[best].region);
    }
    shape.delete();
    ex.delete();
    result.faces = index; result.free = 0;
    for (const list of edgeUse.values()) for (const item of list) if (item.count === 1) {
      result.free++;
      for (const f of item.faces) if (faceRegion[f] !== undefined) result.bad.add(faceRegion[f]);
    }
    if (result.free) result.valid = false;
    // Área por región (las dos mitades de una esfera o toro completos comparten región).
    let total = 0; for (const a of regionArea) total += a;
    const slack = 5e-4 * total + 0.02;
    const byRegion = new Map();
    brep.faces.forEach((f, i) => { if (!f.loops.length) return; const e = byRegion.get(f.region) || { area: 0, seen: true }; e.area += got[i]; e.seen = e.seen && !!seen[i]; byRegion.set(f.region, e); });
    for (const [region, e] of byRegion) {
      const want = regionArea[region];
      if (e.seen && (e.area < 0 || e.area > 2.5 * want + slack)) (result.gross = result.gross || new Set()).add(region);
      if (!e.seen || Math.abs(e.area - want) > 0.2 * want + slack) { result.bad.add(region); (result.why = result.why || []).push([region, e.seen, +want.toFixed(4), +e.area.toFixed(4)]); }
    }
    result.unmatched = unmatched.length;
    return result;
  } catch (error) {
    result.error = typeof error === 'number' ? 'excepción del motor CAD' : String(error.message || error);
    return result;
  } finally {
    reader.delete(); progress.delete();
    try { oc.FS.unlink('/check.step'); } catch { /* ya eliminado */ }
  }
}
