from Dat import *
from General.Numerical import *
from numpy.core.numeric import NaN

import General.SimpPlot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.12, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(5,4)


def getFirstN(d, iVal):
    for N in sorted(d[iVal]):
        if len((d[iVal][N][INDEX_POLE])) > 0:
            return [float(N)]
    return [0.0]

def plotFirstN(df, index, desc, showMarkers=True):
    ret = plotData("./Data_"+str(df)+"/", index, getFirstN)
    sp.plotSingle("First N with Pole with comparison Threshold df="+str(df*100)+"%", ret[0], ret[1], desc+" Energy Range (Rydbergs)", "N", path="Results/FirstNPole_"+str(df)+".png",markerSz=8 if showMarkers else None,markWithLine=True) 


def getNumPoles(d, iVal):
    return [float(len(d[iVal][65][INDEX_POLE]))]

def plotNumPoles(df, index, desc, showMarkers=True):
    ret = plotData("./Data_"+str(df)+"/", index, getNumPoles)
    sp.plotSingle("Num of Poles found with comparison Threshold df="+str(df*100)+"%", ret[0], ret[1], desc+" Energy Range (Rydbergs)", "Num Poles Found", path="Results/Num_"+str(df)+".png",markerSz=8 if showMarkers else None, markWithLine=True) 

 
def getNumPolesNotLost(d, iVal):
    num = 0.0
    for pole in d[iVal][65][INDEX_POLE]:
        if not pole.isLost:
            num += 1.0
    return [num]

def plotNumPolesNotLost(df, index, desc, showMarkers=True):
    ret = plotData("./Data_"+str(df)+"/", index, getNumPolesNotLost)
    sp.plotSingle("Num of Poles found with comparison Threshold df="+str(df*100)+"%", ret[0], ret[1], desc+" Energy Range (Rydbergs)", "Num Poles Found", path="Results/Num_"+str(df)+".png",markerSz=8 if showMarkers else None, markWithLine=True) 


#Assume the first pole is the actual pole
lastPole = None
def getValue(d, iVal):
    global lastPole
    if len((d[iVal][65][INDEX_POLE])) > 0:
        c = Compare(0.001)
        for pole in d[iVal][65][INDEX_POLE]:
            if lastPole is None or c.complexCompare(lastPole.E, pole.E):
                val = pole.E
                lastPole = pole
                return [val.real, val.imag]
    return [NaN, NaN]

def plotValues(df, index, desc, imag, showMarkers=True):
    ret = plotData("./Data_"+str(df)+"/", index, getValue)
    sp.plotSingle(str("Imag" if imag else "Real") + " values of Poles found with comp thres df="+str(df*100)+"%", ret[0], [ret[1][1]] if imag else [ret[1][0]], desc+" Energy Range (Rydbergs)", "Energy (Rydbergs)", path="Results/Values_"+str("Imag_" if imag else "Real_")+str(df)+".png",markerSz=8 if showMarkers else None,markWithLine=True) 


def plotData(dataPath, pathIndex, fun):
    d = getData(dataPath,pathIndex).getData()
    xs = []
    ys_s = []
    first = True
    for iVal in sorted(d.keys()):
      try:
          ys = fun(d,iVal)
          if first:
            for i in range(len(ys)): 
              ys_s.append([])
          for i in range(len(ys_s)):
            ys_s[i].append(ys[i])
          first = False
      except:
          pass #Most likely an issue with the calculation file (eg Too many iteration in Lagurre. Only want to add x if not the case.
      else:
          xs.append(float(iVal)*STEP_TO_RYD)
    return(xs, ys_s)
