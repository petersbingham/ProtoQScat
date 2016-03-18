from runArgsPlot import *
from runBase import *
import Scattering.Stran as S

anaSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]
fitSignList = [[1.0,1.0]]
ratSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]

sm.setSubPlots(len(anaSignList),len(fitSignList)*len(ratSignList)+1,getParameterDesc(args)+"Cross Sections", "Total Energy (hartrees)", "Cross Section (bohr^2)")

def plotAnaXS(anakCal, anaSmat, signString):
    anaXSmat = S.XSmat(S.Tmat(anaSmat), anakCal)
    doMatPlot(args, anaXSmat, False, "Analytical", signString)

def plotRatXS(ratkCal, ratSmat, signString):
    ratXSmat = S.XSmat(S.Tmat(ratSmat), ratkCal)
    doMatPlot(args, ratXSmat, False, "Approximated", signString)
    
dokSignIt(args, anaSignList, fitSignList, ratSignList, plotAnaXS, plotRatXS)