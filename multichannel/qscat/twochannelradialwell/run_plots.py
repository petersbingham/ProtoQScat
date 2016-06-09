from runargsplot import *

anaSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]
fitSignList = [[1.0,1.0]]
ratSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]
IMAG = False

sm.setExtents([1,8],[-2,2])
sm.setParams(0.09, 0.07, 0.99, 0.92, 0.05, 0.2)
sm.setSize(18.5, 10.5)
sm.setSubPlots(len(anaSignList),len(fitSignList)*len(ratSignList)+1,getParameterDesc(args)+"S matrices", "Total Energy (hartrees)", "Imag" if IMAG else "Real")

def plotAnaS(anakCal, anaSmat, signString):
    doMatPlot(args, anaSmat, IMAG, "Analytical,  Imag Offset:" + str(args.eneComplex_), signString)

def plotRatS(ratkCal, ratSmat, signString):
    doMatPlot(args, ratSmat, IMAG, "Approximated,  Imag Offset:" + str(args.eneComplex_), signString)
    
dokSignIt(args, anaSignList, fitSignList, ratSignList, plotAnaS, plotRatS)
