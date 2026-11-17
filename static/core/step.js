// Escritura de STEP (AP214) con caras analíticas a partir del B-rep.
import { perp, dot, sub, mul, add, cross, unit } from './fit.js';

function num(x) {
  if (!Number.isFinite(x)) throw new Error('Número no válido en la geometría.');
  let s = String(+x.toPrecision(15));
  if (s.includes('e')) { let [m, e] = s.split('e'); if (!m.includes('.')) m += '.'; return m + 'E' + e; }
  return s.includes('.') ? s : s + '.';
}

export function writeStep(brep, name = 'pieza', uncertainty = 1e-6, options = {}) {
  const lines = []; let id = 0;
  const add_ = text => { lines.push(`#${++id}=${text};`); return id; };
  const pt = p => add_(`CARTESIAN_POINT('',(${num(p[0])},${num(p[1])},${num(p[2])}))`);
  const dir = d => add_(`DIRECTION('',(${num(d[0])},${num(d[1])},${num(d[2])}))`);
  const axis2 = (o, z, x) => add_(`AXIS2_PLACEMENT_3D('',#${pt(o)},#${dir(z)},#${dir(x)})`);
  const safe = String(name).replace(/[^A-Za-z0-9_ .-]/g, '_');

  const appCtx = add_(`APPLICATION_CONTEXT('core data for automotive mechanical design processes')`);
  add_(`APPLICATION_PROTOCOL_DEFINITION('international standard','automotive_design',2000,#${appCtx})`);
  const prodCtx = add_(`PRODUCT_CONTEXT('',#${appCtx},'mechanical')`);
  const product = add_(`PRODUCT('${safe}','${safe}','',(#${prodCtx}))`);
  add_(`PRODUCT_RELATED_PRODUCT_CATEGORY('part',$,(#${product}))`);
  const formation = add_(`PRODUCT_DEFINITION_FORMATION('','',#${product})`);
  const defCtx = add_(`PRODUCT_DEFINITION_CONTEXT('part definition',#${appCtx},'design')`);
  const definition = add_(`PRODUCT_DEFINITION('design','',#${formation},#${defCtx})`);
  const shape = add_(`PRODUCT_DEFINITION_SHAPE('','',#${definition})`);
  const mm = add_(`(LENGTH_UNIT()NAMED_UNIT(*)SI_UNIT(.MILLI.,.METRE.))`);
  const rad = add_(`(NAMED_UNIT(*)PLANE_ANGLE_UNIT()SI_UNIT($,.RADIAN.))`);
  const sr = add_(`(NAMED_UNIT(*)SI_UNIT($,.STERADIAN.)SOLID_ANGLE_UNIT())`);
  const unc = add_(`UNCERTAINTY_MEASURE_WITH_UNIT(LENGTH_MEASURE(${num(uncertainty)}),#${mm},'distance_accuracy_value','confusion accuracy')`);
  const ctx = add_(`(GEOMETRIC_REPRESENTATION_CONTEXT(3)GLOBAL_UNCERTAINTY_ASSIGNED_CONTEXT((#${unc}))GLOBAL_UNIT_ASSIGNED_CONTEXT((#${mm},#${rad},#${sr}))REPRESENTATION_CONTEXT('Context #1','3D Context with UNIT and UNCERTAINTY'))`);
  const origin = axis2([0, 0, 0], [0, 0, 1], [1, 0, 0]);

  const vertexIds = brep.vertices.map(p => add_(`VERTEX_POINT('',#${pt(p)})`));
  const edgeIds = brep.edges.map(e => {
    const c = e.curve; let curve;
    if (c.type === 'line') curve = add_(`LINE('',#${pt(c.p)},#${add_(`VECTOR('',#${dir(c.d)},1.)`)})`);
    else if (c.type === 'circle') curve = add_(`CIRCLE('',#${axis2(c.c, c.a, c.x)},${num(c.r)})`);
    else if (c.type === 'bspline') {
      const pts = c.ctrl.map(p => '#' + pt(p)).join(',');
      curve = add_(`B_SPLINE_CURVE_WITH_KNOTS('',${c.degree},(${pts}),.UNSPECIFIED.,.F.,.F.,(${c.mults.join(',')}),(${c.knots.map(num).join(',')}),.UNSPECIFIED.)`);
    } else throw new Error('Curva no admitida.');
    return add_(`EDGE_CURVE('',#${vertexIds[e.v0]},#${vertexIds[e.v1]},#${curve},.T.)`);
  });

  const surface = (S, f) => {
    switch (S.type) {
      case 'plane': { const o = mul(S.n, S.d); return add_(`PLANE('',#${axis2(o, S.n, perp(S.n))})`); }
      case 'cylinder': return add_(`CYLINDRICAL_SURFACE('',#${axis2(S.p, S.a, perp(S.a))},${num(S.r)})`);
      case 'sphere': {
        // Medias esferas completas (un único círculo máximo): eje Z como en su contorno.
        const loops = f.loops, e0 = loops.length === 1 && loops[0].length === 1 ? brep.edges[loops[0][0].edge] : null;
        if (e0 && e0.curve.type === 'circle' && Math.abs(e0.curve.r - S.r) < 1e-9) return add_(`SPHERICAL_SURFACE('',#${axis2(S.c, [0, 0, 1], [1, 0, 0])},${num(S.r)})`);
        // Resto: los polos de la parametrización, lejos de la cara (cara centrada en el ecuador).
        let m = [0, 0, 0];
        for (const loop of loops) for (const ref of loop) for (const v of [brep.edges[ref.edge].v0, brep.edges[ref.edge].v1]) { const u = unit(sub(brep.vertices[v], S.c)); m = [m[0] + u[0], m[1] + u[1], m[2] + u[2]]; }
        m = unit(m);
        return add_(`SPHERICAL_SURFACE('',#${axis2(S.c, perp(m), m)},${num(S.r)})`);
      }
      case 'torus': return add_(`TOROIDAL_SURFACE('',#${axis2(S.c, S.a, perp(S.a))},${num(S.R)},${num(S.r)})`);
      case 'extrusion': {
        const c = S.curve, pts = c.ctrl.map(p => '#' + pt(p)).join(',');
        const curve = add_(`B_SPLINE_CURVE_WITH_KNOTS('',${c.degree},(${pts}),.UNSPECIFIED.,.F.,.F.,(${c.mults.join(',')}),(${c.knots.map(num).join(',')}),.UNSPECIFIED.)`);
        return add_(`SURFACE_OF_LINEAR_EXTRUSION('',#${curve},#${add_(`VECTOR('',#${dir(S.d)},1.)`)})`);
      }
      case 'cone': {
        // Situar el origen en el extremo ancho de la cara para que el radio sea positivo.
        let h = 0;
        for (const loop of f.loops) for (const ref of loop) for (const v of [brep.edges[ref.edge].v0, brep.edges[ref.edge].v1]) h = Math.max(h, dot(sub(brep.vertices[v], S.apex), S.a));
        if (!(h > 0)) h = 1;
        return add_(`CONICAL_SURFACE('',#${axis2(add(S.apex, mul(S.a, h)), S.a, perp(S.a))},${num(h * Math.tan(S.alpha))},${num(S.alpha)})`);
      }
    }
    throw new Error('Superficie no admitida.');
  };

  const faceIds = brep.faces.map(f => {
    if (!f.loops.length) return 0;
    const bounds = f.loops.map((loop, i) => {
      const oriented = loop.map(r => '#' + add_(`ORIENTED_EDGE('',*,*,#${edgeIds[r.edge]},${r.forward ? '.T.' : '.F.'})`)).join(',');
      const el = add_(`EDGE_LOOP('',(${oriented}))`);
      return '#' + add_(`${i === (f.outer || 0) ? 'FACE_OUTER_BOUND' : 'FACE_BOUND'}('',#${el},.T.)`);
    });
    return add_(`ADVANCED_FACE('',(${bounds.join(',')}),#${surface(f.S, f)},${f.sign > 0 ? '.T.' : '.F.'})`);
  });

  const solids = [], shellIds = brep.shells.map(shell => (shell.faces.length ? add_(`CLOSED_SHELL('',(${shell.faces.map(i => '#' + faceIds[i]).join(',')}))`) : 0));
  const voidFlag = options.voidFlag ?? '.F.';
  brep.shells.forEach((shell, c) => {
    if (!shellIds[c] || shell.parent >= 0) return;
    const voids = brep.shells.map((x, i) => (x.parent === c && shellIds[i] ? i : -1)).filter(i => i >= 0);
    if (!voids.length) { solids.push(add_(`MANIFOLD_SOLID_BREP('',#${shellIds[c]})`)); return; }
    const oriented = voids.map(i => '#' + add_(`ORIENTED_CLOSED_SHELL('',*,#${shellIds[i]},${voidFlag})`)).join(',');
    solids.push(add_(`BREP_WITH_VOIDS('',#${shellIds[c]},(${oriented}))`));
  });
  const rep = add_(`ADVANCED_BREP_SHAPE_REPRESENTATION('',(#${origin},${solids.map(s => '#' + s).join(',')}),#${ctx})`);
  add_(`SHAPE_DEFINITION_REPRESENTATION(#${shape},#${rep})`);

  const stamp = new Date().toISOString().slice(0, 19);
  return `ISO-10303-21;\nHEADER;\nFILE_DESCRIPTION(('STL convertido a STEP con superficies analiticas'),'2;1');\nFILE_NAME('${safe}.step','${stamp}',(''),(''),'Conversor 3D','Conversor 3D','');\nFILE_SCHEMA(('AUTOMOTIVE_DESIGN { 1 0 10303 214 1 1 1 1 }'));\nENDSEC;\nDATA;\n${lines.join('\n')}\nENDSEC;\nEND-ISO-10303-21;\n`;
}
void cross; void unit;
