import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../Utilities/General')
import SimpPlot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.16, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(6,3)
import numpy as np

if sys.argv[1] == 'both':
  real = True
  imag = True
  axlabel = "Domain"
elif sys.argv[1] == 'real':
  real = True
  imag = False
  axlabel = "Range real"
elif sys.argv[1] == 'imag':
  real = False
  imag = True
  axlabel = "Range imag"
  
val = float(sys.argv[2])

vals = np.arange(-5, 5, 0.05)
  
fr = lambda x:x.real
fi = lambda x:x.imag

def _val(X,Y,Th):
  return np.sqrt((X+Y*1.0j)-Th)
  
desc = [0.0,3.0]
Zposs = []
for Th in desc:
  for sign in [-1.0, 1.0]:
    if real:
      Zposs.append(map(fr, sign*_val(val,vals,Th)))
    if imag:
      Zposs.append(map(fi, sign*_val(val,vals,Th)))

sp.plotSingle("Cross Section",vals,Zposs,xlabel="real",ylabel="",drawAxes=True)