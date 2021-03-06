from runargsplot import *
from runbase import *
from scattering.stran import *

anaSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]
fitSignList = [[1.0,1.0]]
ratSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]

sm.setExtents([1,8],[0.005,2500])
sm.setParams(0.09, 0.07, 0.99, 0.92, 0.05, 0.2)
sm.setSize(18.5, 10.5)
sm.setSubPlots(len(anaSignList),len(fitSignList)*len(ratSignList)+1,getParameterDesc(args)+"Cross Sections", "Total Energy (hartrees)", "Cross Section (bohr^2)")

def plotAnaXS(anakCal, anaSmat, signString):
    anaXSmat = T_to_XS(S_to_T(anaSmat), anakCal)
    doMatPlot(args, anaXSmat, False, "Analytical", signString)

def plotRatXS(ratkCal, ratSmat, signString):
    ratXSmat = T_to_XS(S_to_T(ratSmat), ratkCal)
    doMatPlot(args, ratXSmat, False, "Approximated", signString)
    
dokSignIt(args, anaSignList, fitSignList, ratSignList, plotAnaXS, plotRatXS)