import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../Utilities/General')

import cmath
import simpplot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.16, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(6,3)

def pfun(E):
  return cmath.sqrt(E-2)

E_REALS = []
P_REALS_POS = []
P_IMAGS_POS = []
P_REALS_NEG = []
P_IMAGS_NEG = []

for i in range(0,500,1):
  E = (float(i)*4.0)/500.0
  p = pfun(E)
  E_REALS.append(E)
  P_REALS_POS.append(p.real)
  P_IMAGS_POS.append(p.imag)
  P_REALS_NEG.append(-p.real)
  P_IMAGS_NEG.append(-p.imag)

sp.plotSingle("p plane",E_REALS,[P_REALS_POS,P_REALS_NEG],"Real Energy","Real Momentum",path="Results/real.png",legends=["+","-"])
sp.plotSingle("p plane",E_REALS,[P_IMAGS_POS,P_IMAGS_NEG],"Real Energy","Imag Momentum",path="Results/imag.png",legends=["+","-"])
