import sys
import os
import argparse

parentArgs = argparse.ArgumentParser(description="Single Channel Radial Well. Data Analysis.")
parentArgs.add_argument("a_", help="Width of Well.", type=float)
parentArgs.add_argument("results_", help="1-Root Table, 2-Pole Table, 3-Root Total Plot, 4-Pole Total Plot, 5-Root Scat Plot", type=int)
parentArgs.add_argument("range_", help="1-Vary kmin, 2-Vary kmax", type=int)
parentArgs.add_argument("lockVal_", help="Value of non-varying parameter.", type=float)
parentArgs.add_argument("chopEnd_", help="For scattering number of spurious states to remove.", nargs='?', type=int, default=None)
parentArgs.add_argument("logx_", help="Plot logarithmic x axis.", nargs='?', type=bool, default=False)
args = parentArgs.parse_args()

base_kvar =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base_kvar+'./../..')
from Base import *

import General.SimpPlot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.12, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(5,4)
sp.turnOffColourCycle()

l = 0
if args.range_==2:
  kmin = args.lockVal_
else:
  kmax = args.lockVal_
df = 0.001
N = 65
astr = "%.1f"%args.a_
resultsPath = base_kvar+"/a=" + astr + "/Results/"

p = Printer(True, False)
f = Filter()

d = getData(base_kvar+'/../../../Data/'+'0.001/'+astr+'/')
if args.range_==1:
  kmin_s = None
  kmax_s = [kmax]  
else: 
  kmin_s = [kmin]
  kmax_s = None
f.selectWhereIn(d.getData(),a_s=[args.a_],l_s=[l],kmin_s=kmin_s,kmax_s=kmax_s,df_s=[df],N_s=[N])
if args.range_==1:
  plotLvl = LVL_KMIN
  f.selectWhereLessThan(d.getData(),kmin=kmax)
else:
  plotLvl = LVL_KMAX
  f.selectWhereGreaterThan(d.getData(),kmax=kmin)
  
f.clean(d.getData())
#d.setComparator(num.Compare(0.001))
d.setComparator(num.Compare(0.003))

if args.results_==2 or args.results_==4:
  searchFun = d.searchEForPole
  title_base = "a=" + astr + ", Poles PADE("
else:
  searchFun = d.searchEForRoot
  title_base = "a=" + astr + ", Roots PADE("
if args.range_==1:
  plotLvl = LVL_KMIN
else:
  plotLvl = LVL_KMAX

if args.range_==1:
  range = d.getRange(LVL_KMIN)
  title = title_base + "N" + str(N) + "_kmin" + str(range[0]) + "-" + str(range[1]) + "_kmax" + str(kmax) + ")"
else:
  range = d.getRange(LVL_KMAX)
  title = title_base + "N" + str(N) + "_kmin" + str(kmin) + "_kmax" + str(range[0]) + "-" + str(range[1]) + ")"
  
foundData = d.dictSearch(searchFun, searchDict[args.a_], LVL_V)
if args.results_==1 or args.results_==2:
  sys.stdout = open(resultsPath+"Table "+title+".txt", 'w')
  p.tabulateListSearch_E(foundData,LVL_V)
else:
  if args.logx_==1:
    xLabel = "kmin"
  else:
    xLabel = "kmax"
  if args.results_==3 or args.results_==4:
    yLabel = "Num 1st Bound States found"
    xs, ys = p.totalListSearch_E(foundData,LVL_V,plotLvl)
    sp.plotSingle(title, xs, [ys], xLabel, yLabel, logx=args.logx_, path=resultsPath+"Chart Total "+title+".png")
  elif args.results_==5:
    yLabel = "1st Bound State Energy (hartrees)"
    xs, ys = p.plotListSearch_E(foundData,LVL_V,plotLvl)
    sp.setSubPlotParameters(left=0.15)
    if args.chopEnd_ is not None:
      ys = ys[:-args.chopEnd_]
    sp.plotSingle(title, xs, ys, xLabel, yLabel, logx=args.logx_, markerSz=5,path=resultsPath+"Chart Scat "+title+".png")  

  
  
  
  
