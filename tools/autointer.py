"""Lista las parejas de caras que se cruzan (autointersecciones) en un STEP: python tools/autointer.py pieza.step [n]"""
import sys
from OCP.STEPControl import STEPControl_Reader
from OCP.BRepAlgoAPI import BRepAlgoAPI_Check
from OCP.Message import Message, Message_Gravity, Message_ProgressRange
from OCP.TopAbs import TopAbs_FACE, TopAbs_EDGE, TopAbs_VERTEX
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
for pr in Message.DefaultMessenger_s().Printers(): pr.SetTraceLevel(Message_Gravity.Message_Fail)
names = {0: 'plano', 1: 'cilindro', 2: 'cono', 3: 'esfera', 4: 'toro'}
def desc(sh):
    t = sh.ShapeType()
    if t == TopAbs_FACE:
        f = TopoDS.Face(sh); g = GProp_GProps(); BRepGProp.SurfaceProperties_s(f, g); c = g.CentreOfMass()
        return f"{names.get(int(BRepAdaptor_Surface(f).GetType()), 'otra')}(área {g.Mass():.3f} en {c.X():.2f},{c.Y():.2f},{c.Z():.2f})"
    if t == TopAbs_VERTEX:
        from OCP.BRep import BRep_Tool
        v = TopoDS.Vertex(sh); q = BRep_Tool.Pnt_s(v); return f"vértice({q.X():.3f},{q.Y():.3f},{q.Z():.3f} tol {BRep_Tool.Tolerance_s(v):.4f})"
    g = GProp_GProps(); BRepGProp.LinearProperties_s(sh, g); c = g.CentreOfMass()
    return f"{str(t).split('_')[-1].lower()}(long {g.Mass():.4f} en {c.X():.2f},{c.Y():.2f},{c.Z():.2f})"
r = STEPControl_Reader(); r.ReadFile(sys.argv[1]); r.TransferRoots(); s = r.OneShape()
c = BRepAlgoAPI_Check(s, True, True, Message_ProgressRange())
n = 0
for res in c.Result():
    if (sys.argv[3] if len(sys.argv) > 3 else 'SelfIntersect') not in str(res.GetCheckStatus()): continue
    n += 1
    if n <= (int(sys.argv[2]) if len(sys.argv) > 2 else 8):
        print('  ', ' x '.join(desc(x) for x in res.GetFaultyShapes1()))
print('autointersecciones:', n)
