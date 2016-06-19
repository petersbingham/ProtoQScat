import sys
sys.path.append("../../utilities")
import general.simpplot as sp

from math import *

resonances = []

emin = float(sys.argv[1])
emax = float(sys.argv[2])
steps = int(sys.argv[3])

i = 4
while True:
    if i < len(sys.argv)-1:
        res = (float(sys.argv[i]),-2.0*float(sys.argv[i+1]))
        print str(res)
        resonances.append(res)
        i += 2
    else:
        break

def phaseFun(pos, width, ene):
    x = (width/2.0) / (pos-ene)
    return atan(x)

de = (emax-emin) / steps

enes = []
phases = []

ene = emin
for  i in range(steps):
    phase = 0.0
    for res in resonances:
        phase += phaseFun(res[0],res[1],ene)
    enes.append(ene)
    phases.append(phase)
    ene += de

sp.setExtents((0.0,1.6), (-3.5,3.5))
sp.plotSingle("Resonant Phase Shift Plot - " + str(len(resonances))  + " poles",enes,[phases],"Energy (ryds)", "Phase Shift")