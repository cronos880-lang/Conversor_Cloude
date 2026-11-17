"""Caras diminutas y aristas muy cortas (lo que suele romper Parasolid/Onshape): python tools/smallstuff.py a.step"""
import sys, json
from OCP.STEPControl import STEPControl_Reader
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.Message import Message, Message_Gravity
for p in Message.DefaultMessenger_s().Printers(): p.SetTraceLevel(Message_Gravity.Message_Fail)
r = STEPControl_Reader(); r.ReadFile(sys.argv[1]); r.TransferRoots(); s = r.OneShape()
def props(kind, fn):
    e = TopExp_Explorer(s, kind); out = []
    while e.More():
        g = GProp_GProps(); fn(e.Current(), g); out.append(g.Mass()); e.Next()
    return out
fa = props(TopAbs_FACE, BRepGProp.SurfaceProperties_s); ed = props(TopAbs_EDGE, BRepGProp.LinearProperties_s)
print(json.dumps({'caras': len(fa), 'caras<1e-3mm2': sum(a < 1e-3 for a in fa), 'caras<1e-5mm2': sum(a < 1e-5 for a in fa),
  'aristas': len(ed) // 2, 'aristas<0.01mm': sum(l < 0.01 for l in ed) // 2, 'aristas<0.001mm': sum(l < 0.001 for l in ed) // 2, 'arista_min': min(ed)}))
