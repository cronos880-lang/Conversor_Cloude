"""Diagnóstico de caras no válidas de un STEP: python tools/badfaces.py a.step"""
import sys
from OCP.STEPControl import STEPControl_Reader
from OCP.BRepCheck import BRepCheck_Analyzer, BRepCheck_Face, BRepCheck_Wire, BRepCheck_Edge
from OCP.TopoDS import TopoDS_Edge
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_WIRE, TopAbs_EDGE, TopAbs_VERTEX, TopAbs_SHELL
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.BRep import BRep_Tool
from OCP.Message import Message, Message_Gravity
for p in Message.DefaultMessenger_s().Printers(): p.SetTraceLevel(Message_Gravity.Message_Fail)
r = STEPControl_Reader(); r.ReadFile(sys.argv[1]); r.TransferRoots(); shape = r.OneShape()
an = BRepCheck_Analyzer(shape)
names = {0: 'plane', 1: 'cylinder', 2: 'cone', 3: 'sphere', 4: 'torus'}
def statuses(s):
    return [] if BRepCheck_Analyzer(s).IsValid() else ['MAL']
def sub(s, kind):
    e = TopExp_Explorer(s, kind); out = []
    while e.More(): out.append(e.Current()); e.Next()
    return out
n = 0
for f in sub(shape, TopAbs_FACE):
    if BRepCheck_Analyzer(f).IsValid(): continue
    n += 1
    if n > int(sys.argv[2]) if len(sys.argv) > 2 else n > 12: continue
    face = TopoDS.Face(f); t = names.get(int(BRepAdaptor_Surface(face).GetType()), 'other')
    g = GProp_GProps(); BRepGProp.SurfaceProperties_s(f, g); c = g.CentreOfMass()
    issues = set(statuses(f))
    try:
        bf = BRepCheck_Face(face)
        issues.add('F.intersect=' + str(bf.IntersectWires()).split('_')[-1]); issues.add('F.classify=' + str(bf.ClassifyWires()).split('_')[-1]); issues.add('F.orient=' + str(bf.OrientationOfWires()).split('_')[-1])
        for w in sub(f, TopAbs_WIRE):
            bw = BRepCheck_Wire(TopoDS.Wire(w)); e1 = TopoDS_Edge(); e2 = TopoDS_Edge()
            issues.add('W.orient=' + str(bw.Orientation(face)).split('_')[-1]); issues.add('W.closed=' + str(bw.Closed()).split('_')[-1]); issues.add('W.self=' + str(bw.SelfIntersect(face, e1, e2)).split('_')[-1])
    except Exception as ex: issues.add('exc ' + str(ex)[:60])
    issues = set(i for i in issues if not i.endswith('NoError'))
    for w in sub(f, TopAbs_WIRE): issues |= set('W:' + x for x in statuses(w))
    edges = sub(f, TopAbs_EDGE)
    for e in edges: issues |= set('E:' + x for x in statuses(e))
    tolmax = max([BRep_Tool.Tolerance_s(TopoDS.Edge(e)) for e in edges] or [0])
    print(f'{t:9s} area={g.Mass():10.4f} centro=({c.X():.2f},{c.Y():.2f},{c.Z():.2f}) bucles={len(sub(f, TopAbs_WIRE))} aristas={len(edges)} tolmax={tolmax:.2g} {sorted(issues)}')
print('caras malas:', n, 'shell:', [statuses(s) for s in sub(shape, TopAbs_SHELL)][:3])
