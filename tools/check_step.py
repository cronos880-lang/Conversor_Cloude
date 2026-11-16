"""Comprobación independiente de un STEP con Open CASCADE: python tools/check_step.py a.step [volumen] [area]"""
import sys, json, time
from OCP.STEPControl import STEPControl_Reader
from OCP.IFSelect import IFSelect_RetDone
from OCP.BRepCheck import BRepCheck_Analyzer
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE, TopAbs_SOLID, TopAbs_EDGE
from OCP.TopoDS import TopoDS
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.Message import Message, Message_Gravity
for p in Message.DefaultMessenger_s().Printers(): p.SetTraceLevel(Message_Gravity.Message_Fail)
t0 = time.time()
reader = STEPControl_Reader()
out = {'file': sys.argv[1].split('/')[-1]}
if reader.ReadFile(sys.argv[1]) != IFSelect_RetDone:
    out['error'] = 'no se puede leer'; print(json.dumps(out)); sys.exit(1)
reader.TransferRoots(); shape = reader.OneShape()
out['valid'] = bool(BRepCheck_Analyzer(shape).IsValid())
def count(kind):
    e = TopExp_Explorer(shape, kind); n = 0
    while e.More(): n += 1; e.Next()
    return n
out['solids'] = count(TopAbs_SOLID); out['faces'] = count(TopAbs_FACE)
props = GProp_GProps(); BRepGProp.VolumeProperties_s(shape, props); out['volume'] = round(props.Mass(), 3)
props = GProp_GProps(); BRepGProp.SurfaceProperties_s(shape, props); out['area'] = round(props.Mass(), 3)
names = {0: 'plane', 1: 'cylinder', 2: 'cone', 3: 'sphere', 4: 'torus'}
types = {}; e = TopExp_Explorer(shape, TopAbs_FACE); bad = 0
while e.More():
    f = TopoDS.Face(e.Current()); t = int(BRepAdaptor_Surface(f).GetType()); k = names.get(t, 'other'); types[k] = types.get(k, 0) + 1
    if not BRepCheck_Analyzer(f).IsValid(): bad += 1
    e.Next()
out['types'] = types; out['bad_faces'] = bad
if len(sys.argv) > 2:
    v = float(sys.argv[2]); out['vol_err_%'] = round(100 * (out['volume'] - v) / v, 3)
if len(sys.argv) > 3:
    a = float(sys.argv[3]); out['area_err_%'] = round(100 * (out['area'] - a) / a, 3)
out['read_s'] = round(time.time() - t0, 2)
print(json.dumps(out))
