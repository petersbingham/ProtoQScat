from Dat import *
from General.Numerical import *
from numpy.core.numeric import NaN

import General.SimpPlot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.12, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(8,4)


#Ns = [9,17,33,65]

lastPole = None
def getValue(d, kmax, Nrange):
    global lastPole
    ys_real = []
    ys_imag = []
    for N in range(Nrange[0],Nrange[1]+1):
        printed = False
        if len((d[kmax][N][INDEX_POLE])) > 0:
            c = Compare(0.001)
            for pole in d[kmax][N][INDEX_POLE]:
                if lastPole is None or c.complexCompare(lastPole.E, pole.E):
                    lastPole = pole
                    val = pole.E
                    ys_real.append(val.real)
                    ys_imag.append(val.imag)
                    printed = True
                    break
        if not printed:
          print N
    return ys_real, ys_imag

def getDiff(d, kmax, Nrange):
    global lastPole
    ys_real = []
    ys_imag = []
    for N in range(Nrange[0],Nrange[1]+1):
        printed = False
        if len((d[kmax][N][INDEX_POLE])) > 0:
            c = Compare(0.01)
            for pole in d[kmax][N][INDEX_POLE]:
                if lastPole is None or c.complexCompare(lastPole.E, pole.E):
                    lastPole = pole
                    for root in d[kmax][N-1][INDEX_POLE]:
                        if c.complexCompare(root.E, pole.E):                        
                            val = pole.E - root.E
                            ys_real.append(val.real)
                            ys_imag.append(val.imag)
                            printed = True
                            break
                    if printed:
                      break
        if not printed:
          print N
    return ys_real, ys_imag

def plotValues(df, index, kmax, Nrange, real):
    showMarkers = True
    ret = plotData("./Data_suc_"+str(df)+"/", index, getValue, kmax, Nrange)
    if real:
        desStr = "Real"
    else:
        desStr = "Imag"
    sp.plotSingle(desStr + " Value of Poles found with comparison Threshold df="+str(df*100)+"%, Emax="+stepToRydStr(kmax)+"Ryds", ret[0], ret[1] if real else ret[2], "N", desStr+" Energy (Rydbergs)", path="Results/Values_"+desStr+"_"+str(df)+".png",markerSz=8 if showMarkers else None,markWithLine=True,logy=False) 

def plotDiffs(df, index, kmax, Nrange, real):
    showMarkers = True
    ret = plotData("./Data_suc_"+str(df)+"/", index, getDiff, kmax, Nrange)
    if real:
        desStr = "Real"
    else:
        desStr = "Imag"
    sp.plotSingle(desStr + " root differences for Poles found with comparison threshold df="+str(df*100)+"%, Emax="+stepToRydStr(kmax)+"Ryds", ret[0], ret[1] if real else ret[2], "N", desStr+" Energy Diff (Rydbergs)", path="Results/DiffValues_"+desStr+"_"+str(df)+".png",markerSz=8 if showMarkers else None,markWithLine=True,logy=True) 

def plotData(dataPath, pathIndex, fun, kmax, Nrange):
    d = getData(dataPath,pathIndex).getData()
    xs = [float(N) for N in range(Nrange[0],Nrange[1]+1)]
    ys_real_s = []
    ys_imag_s = []
    first = True
    ys_real, ys_imag = fun(d, kmax, Nrange)
    ys_real_s.append(ys_real)
    ys_imag_s.append(ys_imag)
    print "***** " + str(len(xs)) + "   " + str(len(ys_real_s[0]))
    return(xs, ys_real_s, ys_imag_s)
