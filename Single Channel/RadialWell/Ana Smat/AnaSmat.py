import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../Utilities/General')

POINTS = 800000

import cmath
import sympy.mpmath as mpm
import SimpPlot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.16, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(6,3)

from AnaSmatFuns import *

if len(sys.argv) < 5:
  print "Not enough arguments"
else:
  t = int(sys.argv[1])
  a = float(sys.argv[2])
  V = float(sys.argv[3])
  start = complex(sys.argv[4])
  if len(sys.argv) == 6:
    end = float(sys.argv[5])
  
def denumWrap(k):
  return denum(k,a,V)
  
if t == 1:
  print mpm.findroot(denumWrap, start)
else:
  xs = []
  Ss_real = []
  Ss_imag = []
  if t == 2:
    k = start
    dk = (end-start.real) / POINTS
    for i in range(POINTS):
      xs.append(k.real)
      calValues(k, a, V, Ss_real, Ss_imag)
      k += dk
    name = "Analytical S matrix for Radial Well, V0="+str(V)+", a="+str(a)+", k.imag="+str(k.imag)
    sp.plotSingle(name,xs,[Ss_real,Ss_imag],"k.real","",legends=["S.real","S.imag"],path="Results/"+name+".png")
  elif t == 3:
    k = start
    dk = ((end-start.imag) / POINTS) * 1.0j
    for i in range(POINTS):
      xs.append(k.imag)
      calValues(k, a, V, Ss_real, Ss_imag)
      k += dk
    name = "Analytical S matrix for Radial Well, V0="+str(V)+", a="+str(a)+", k.real="+str(k.real)
    sp.plotSingle(name,xs,[Ss_real,Ss_imag],"k.imag","",legends=["S.real","S.imag"],path="Results/"+name+".png")
    
    