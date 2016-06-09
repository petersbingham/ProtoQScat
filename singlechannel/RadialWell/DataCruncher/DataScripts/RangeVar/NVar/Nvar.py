import sys
import os
import argparse

parentArgs = argparse.ArgumentParser(description="Single Channel Radial Well. Data Analysis.")
parentArgs.add_argument("a_", help="Width of Well.", type=float)
parentArgs.add_argument("results_", help="1-Root Table, 2-Pole Table, 3-Root Total Plot, 4-Pole Total Plot, 5-Root Scat Plot, 6-Pole Scat Plot", nargs='?', type=int, default=1)
parentArgs.add_argument("kmin_", help="Value of kmin", nargs='?', type=float, default=0.01)
parentArgs.add_argument("chopEnd_", help="For scattering number of spurious states to remove.", nargs='?', type=int, default=None)
parentArgs.add_argument("logx_", help="Plot logarithmic x axis.", nargs='?', type=bool, default=False)
args = parentArgs.parse_args()

base_kvar =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base_kvar+'./../..')
from Base import *

import general.simpplot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.12, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(5,4)

l = 0
kmin = args.kmin_
df = 0.001

p = Printer(True, False)
f = Filter()
astr = "%.1f"%args.a_
resultsPath = base_kvar+"/a=" + astr + "/Results/"

def _uniformVals(xs_s, ys_s):
  minx = 999
  maxx = 0
  
  newys = []
  for xs in xs_s:
    if len(xs) > 0:
      if xs[0] < minx:
        minx = xs[0]
      i = len(xs)-1
      if xs[i] > maxx:
        maxx = xs[i]
      
  newxs = [minx+0.5*i for i in range(int((maxx-minx)*2+1.1))]
  newys_s = []
  for ys in ys_s:
    newys_s.append([0.0]*len(newxs)) #This needs to be a float or we get buffer issues with the nparray when plotting.
    
  for i in range(len(xs_s)):
    for j in range(len(xs_s[i])):
      nj = newxs.index(xs_s[i][j])
      newys_s[i][nj] = ys_s[i][j]
  
  return newxs, newys_s

ys_s = []
xs_s = []
legends = []
oldxs = None
kRange = None
for N in [3,5,9,17,33,65]:
  d = getData(base_kvar+'/../../../Data/'+'0.001/'+astr+'/')

  f.selectWhereIn(d.getData(),a_s=[args.a_],l_s=[l],kmin_s=[kmin],df_s=[df],N_s=[N])
  f.selectWhereGreaterThan(d.getData(),kmax=kmin)
    
  f.clean(d.getData())
  d.setComparator(num.Compare(0.003))

  if args.results_==1 or args.results_==3 or args.results_==5:
    searchFun = d.searchEForRoot
    title_base = "a=" + astr + ", Roots PADE("
  else:
    searchFun = d.searchEForPole
    title_base = "a=" + astr + ", Poles PADE("
  plotLvl = LVL_KMAX

  newkRange = d.getRange(LVL_KMAX)
  if kRange is None:
    kRange = newkRange
  else:
    if newkRange[0] < kRange[0]:
      kRange[0] = newkRange[0]
    if newkRange[1] > kRange[1]:
      kRange[1] = newkRange[1]
  
  foundData = d.dictSearch(searchFun, searchDict[args.a_], LVL_V)

  xLabel = "kmax"
  
  title_end = "_kmin"+str(kmin)+"_kmax"+str(kRange[0])+"-"+str(kRange[1]) + ")"
  
  if args.results_==1 or args.results_==2:
    title = title_base + "N" + str(N) + title_end
    sys.stdout = open("Results/Table "+title+".txt", 'w')
    p.tabulateListSearch_E(foundData,LVL_V)  
  elif args.results_==3 or args.results_==4:
    yLabel = "Num 1st Bound States found"
    xs, ys = p.totalListSearch_E(foundData,LVL_V,plotLvl)
    legends.append("N="+str(N))
    xs_s.append(xs)
    ys_s.append(ys)
  else:
    title = title_base + "N" + str(N) + title_end
    yLabel = "1st Bound State Energy (hartrees)"
    xs, ys = p.plotListSearch_E(foundData,LVL_V,plotLvl)
    sp.setSubPlotParameters(left=0.18)
    if args.chopEnd_ is not None:
      ys = ys[:-args.chopEnd_]
    sp.turnOffColourCycle()
    sp.plotSingle(title, xs, ys, xLabel, yLabel, logx=args.logx_, markerSz=5,path="Results/Chart Scat "+title+".png")  

if args.results_==3 or args.results_==4:
  title = title_base + "N_3-65" + title_end
  xs, ys_s = _uniformVals(xs_s, ys_s)
  sp.plotSingle(title, xs, ys_s, xLabel, yLabel, legends=legends, logx=args.logx_, path="Results/Chart Total "+title+".png")
      
    
  
  
  
  
  
