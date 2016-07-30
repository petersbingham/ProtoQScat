import math
import numpy as np
import scipy as sp
from general.qstype import *
 
RANGE_FCT = 1.0 
 
def wid(imag):
    return -imag*2.0
 
def resonant(resData, ene):
    return np.arctan((wid(resData.imag)/2.0)/(resData.real-ene))

def bretWig(coeffs, resData, ene):
    return resonant(resData, ene) + np.polyval(coeffs, ene)
    
def inFitRange(resData, ene):
    rangeMin = resData.real-RANGE_FCT*wid(resData.imag)
    rangeMax = resData.real+RANGE_FCT*wid(resData.imag)
    return ene>rangeMin and ene<rangeMax
    
class EPhaseFitter:
    def __init__(self, ratEPhaseMats, resDatas, polyOrder):
        self.ratEPhaseMats = ratEPhaseMats
        self.resDatas = resDatas
        self.polyOrder = polyOrder
        
    def getCoeffs(self):
        self._extractEnergyRanges()
        return self._doFits()
        
    def _extractEnergyRanges(self):
        self.energyRanges = [[]]*len(self.resDatas)
        for ene in self.ratEPhaseMats:
            for i in range(len(self.resDatas)):
                resData = self.resDatas[i]
                if inFitRange(resData, ene):
                    self.energyRanges[i].append(ene)
        
    def _doFits(self):
        coeffs = []
        for i in range(len(self.energyRanges)):
            startCoeffs = np.array([1.0]*(self.polyOrder+1))
            #startCoeffs = np.array([1000.0,10,-7.0])
            self.currentResData = self.resDatas[i]
            erangeList = self.energyRanges[i]
            eRange = np.array(erangeList)
            epRange = np.array([QStrace(self.ratEPhaseMats[ene]) for ene in erangeList])
            args=(eRange, epRange)
            coeff = self._minimize(startCoeffs, args)
            coeffs.append(coeff)
        return coeffs
    
    def _leastsq(self, startCoeffs, args):
        ret = sp.optimize.leastsq(self._bretWigRes_leastsq, startCoeffs, args, ftol=1.49012e-8, xtol=1.49012e-8)
        print ret
        return ret[0]
    
    def _bretWigRes_leastsq(self, coeffs, ene, ep):
        ret1 = ep - bretWig(coeffs, self.currentResData, ene)
        ret2 = np.sqrt(ret1.real**2 + ret1.imag**2)
        return ret2
    
    def _minimize(self, startCoeffs, args):
        ret = sp.optimize.minimize(self._bretWigRes_minimize, startCoeffs, args, method="Nelder-Mead", options={"maxiter":len(startCoeffs)*2000})
        print ret
        return ret.x
    
    def _bretWigRes_minimize(self, coeffs, ene, ep):
        sum = 0.0
        for i in range(len(ene)):
            sum += (ep[i] - bretWig(coeffs, self.currentResData, ene[i])).real**2.0
        #ret = np.sqrt(sum.real**2 + sum.imag**2)
        return sum
