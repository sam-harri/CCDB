# -f -e -m interppointdatatofld:frompts=chan3D_pts.pts chan3D.xml chan3D_pts.fld
import sys
from NekPy.FieldUtils import *

field = Field(sys.argv, force_output=True, error=True)

InputModule.Create("xml",  field, "chan3D.xml").Run()
ProcessModule.Create("interppointdatatofld", field, frompts="chan3D_pts.pts").Run()
OutputModule.Create("fld", field, "chan3D_pts.fld").Run()
