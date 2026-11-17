"""Auditoría de un STEP con los criterios que suelen hacer fallar la importación en Onshape.
python tools/audit.py pieza.step [pieza.stl]
Mide, además de la validez: aristas fuera de la pieza, holgura arista-cara, aristas y caras diminutas."""
import sys, json, math
from OCP.STEPControl import STEPControl_Reader
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_SOLID
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.BRepAdaptor import BRepAdaptor_Curve, BRepAdaptor_Surface
from OCP.GeomAPI import GeomAPI_ProjectPointOnSurf
from OCP.Bnd import Bnd_Box
from OCP.BRepBndLib import BRepBndLib
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.Message import Message, Message_Gravity
for p in Message.DefaultMessenger_s().Printers(): p.SetTraceLevel(Message_Gravity.Message_Fail)
names = {0: 'plane', 1: 'cylinder', 2: 'cone', 3: 'sphere', 4: 'torus'}
r = STEPControl_Reader(); r.ReadFile(sys.argv[1]); r.TransferRoots(); shape = r.OneShape()
box = Bnd_Box(); BRepBndLib.Add_s(shape, box); mn, mx = box.CornerMin(), box.CornerMax(); x0, y0, z0, x1, y1, z1 = mn.X(), mn.Y(), mn.Z(), mx.X(), mx.Y(), mx.Z()
diag = math.dist((x0, y0, z0), (x1, y1, z1))
out = {'valid': bool(BRepCheck_Analyzer(shape).IsValid())}
e = TopExp_Explorer(shape, TopAbs_SOLID); n = 0
while e.More(): n += 1; e.Next()
out['solids'] = n
gaps, outside, short, worst = [], 0, 0, []
fx = TopExp_Explorer(shape, TopAbs_FACE)
while fx.More():
    face = TopoDS.Face(fx.Current()); surf = BRep_Tool.Surface_s(face); kind = names.get(int(BRepAdaptor_Surface(face).GetType()), 'other')
    ex = TopExp_Explorer(face, TopAbs_EDGE)
    while ex.More():
        edge = TopoDS.Edge(ex.Current()); ex.Next()
        if BRep_Tool.Degenerated_s(edge): continue
        c = BRepAdaptor_Curve(edge); u0, u1 = c.FirstParameter(), c.LastParameter()
        g = GProp_GProps(); BRepGProp.LinearProperties_s(edge, g); L = g.Mass()
        if L < 1e-3: short += 1
        pts = [c.Value(u0 + (u1 - u0) * k / 8) for k in range(9)]
        if any(p.X() < x0 - 0.01 * diag or p.X() > x1 + 0.01 * diag or p.Y() < y0 - 0.01 * diag or p.Y() > y1 + 0.01 * diag or p.Z() < z0 - 0.01 * diag or p.Z() > z1 + 0.01 * diag for p in pts): outside += 1
        gap = 0
        for p in pts:
            pr = GeomAPI_ProjectPointOnSurf(p, surf)
            if pr.NbPoints(): gap = max(gap, pr.LowerDistance())
        gaps.append(gap)
        if gap > 0.01: worst.append((round(gap, 3), round(L, 3), kind, c.GetType().name if hasattr(c.GetType(), 'name') else str(c.GetType()), [round(v, 2) for v in (pts[4].X(), pts[4].Y(), pts[4].Z())]))
    fx.Next()
faces = []; e = TopExp_Explorer(shape, TopAbs_FACE)
while e.More():
    g = GProp_GProps(); BRepGProp.SurfaceProperties_s(e.Current(), g); faces.append(g.Mass()); e.Next()
gaps.sort()
out.update({'caras': len(faces), 'caras<0.001mm2': sum(a < 1e-3 for a in faces), 'aristas': len(gaps), 'aristas_cortas<0.001': short // 2,
            'aristas_fuera_de_la_pieza': outside // 2, 'holgura_max': round(gaps[-1], 4) if gaps else 0,
            'holgura_p99': round(gaps[int(0.99 * (len(gaps) - 1))], 4) if gaps else 0, 'aristas_holgura>0.01': sum(g > 0.01 for g in gaps), 'aristas_holgura>0.001': sum(g > 0.001 for g in gaps)})
print(json.dumps(out))
for w in sorted(worst, reverse=True)[:int(sys.argv[2]) if len(sys.argv) > 2 else 5]: print('   holgura', w)
