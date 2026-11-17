"""Desviación entre un STEP y el STL original: python tools/deviation.py pieza.step pieza.stl
Tesela cada cara del STEP y mide la distancia de sus puntos a la malla original."""
import sys, json, numpy as np, trimesh
from OCP.STEPControl import STEPControl_Reader
from OCP.BRepMesh import BRepMesh_IncrementalMesh
from OCP.TopExp import TopExp_Explorer
from OCP.TopAbs import TopAbs_FACE
from OCP.TopoDS import TopoDS
from OCP.BRep import BRep_Tool
from OCP.TopLoc import TopLoc_Location
from OCP.BRepAdaptor import BRepAdaptor_Surface
from OCP.BRepGProp import BRepGProp
from OCP.GProp import GProp_GProps
r = STEPControl_Reader(); r.ReadFile(sys.argv[1]); r.TransferRoots(); shape = r.OneShape()
mesh = trimesh.load(sys.argv[2], force='mesh')
BRepMesh_IncrementalMesh(shape, 0.02, False, 0.2, True)
names = {0: 'plane', 1: 'cylinder', 2: 'cone', 3: 'sphere', 4: 'torus'}
e = TopExp_Explorer(shape, TopAbs_FACE); rows = []
while e.More():
    f = TopoDS.Face(e.Current()); loc = TopLoc_Location(); tri = BRep_Tool.Triangulation_s(f, loc)
    if tri is not None and tri.NbNodes() > 0:
        pts = np.array([[p.X(), p.Y(), p.Z()] for p in (tri.Node(i + 1).Transformed(loc.Transformation()) for i in range(tri.NbNodes()))])
        d = trimesh.proximity.closest_point(mesh, pts)[1]
        g = GProp_GProps(); BRepGProp.SurfaceProperties_s(f, g)
        rows.append((float(d.max()), names.get(int(BRepAdaptor_Surface(f).GetType()), 'other'), g.Mass(), pts[d.argmax()].round(2).tolist()))
    e.Next()
rows.sort(reverse=True)
print(json.dumps({'faces': len(rows), 'max_dev': round(rows[0][0], 3), 'faces_over_0.1mm': sum(1 for x in rows if x[0] > 0.1), 'size': mesh.extents.round(1).tolist()}))
for d, t, a, p in rows[:8]: print(f'  {d:8.3f} mm  {t:9s} area={a:9.3f}  en {p}')
