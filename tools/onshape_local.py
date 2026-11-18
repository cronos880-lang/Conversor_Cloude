"""Comprobación local que imita el chequeo de sólidos de Onshape (Parasolid), sin gastar cuota.
python tools/onshape_local.py pieza.step [...]
Usa BRepAlgoAPI_Check de Open CASCADE (autointersecciones y aristas pequeñas)."""
import sys, json, time
from OCP.STEPControl import STEPControl_Reader
from OCP.BRepAlgoAPI import BRepAlgoAPI_Check
from OCP.Message import Message, Message_Gravity, Message_ProgressRange
for pr in Message.DefaultMessenger_s().Printers(): pr.SetTraceLevel(Message_Gravity.Message_Fail)
for f in sys.argv[1:]:
    t0 = time.time(); r = STEPControl_Reader(); r.ReadFile(f); r.TransferRoots(); s = r.OneShape()
    c = BRepAlgoAPI_Check(s, True, True, Message_ProgressRange())
    estados = {}
    try:
        for res in c.Result():
            k = str(res.GetCheckStatus()).split('.')[-1]; estados[k] = estados.get(k, 0) + 1
    except Exception as e: estados['?'] = str(e)[:60]
    print(json.dumps({'pieza': f.split('/')[-1].rsplit('.', 1)[0], 'ok': bool(c.IsValid()), 'faltas': estados, 'segundos': round(time.time() - t0, 1)}), flush=True)
