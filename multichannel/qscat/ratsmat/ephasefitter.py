import math
import numpy as np
import scipy as sp
import general.type_wrap as tw
import time
 
RANGE_FCT = 1.0 
 
def wid(imag):
    return -imag*2.0
 
def resonant(resData, ene):
    return np.arctan((wid(resData.imag)/2.0)/(resData.real-ene))

def bretWig(resData, coeffs, ene):
    return resonant(resData, ene) + np.polyval(coeffs, ene)
    
def inFitRange(resData, ene):
    rangeMin = resData.real-RANGE_FCT*wid(resData.imag)
    rangeMax = resData.real+RANGE_FCT*wid(resData.imag)
    return ene>rangeMin and ene<rangeMax
    
class EPhaseFitter:
    def __init__(self, ratEPhaseMats, resDatas, polyOrder, lockBretWig=True):
        self.ratEPhaseMats = ratEPhaseMats
        self.resDatas = resDatas
        self.polyOrder = polyOrder
        self.lockBretWig = lockBretWig
        
    def getCoeffs(self):
        self._extractEnergyRanges()
        return self._doFits1()
        
    def _extractEnergyRanges(self):
        self.energyRanges = [[]]*len(self.resDatas)
        for ene in self.ratEPhaseMats:
            for i in range(len(self.resDatas)):
                resData = self.resDatas[i]
                if inFitRange(resData, ene):
                    self.energyRanges[i].append(ene)
        
    def _doFits1(self):
        coeffs = []
        for i in range(len(self.energyRanges)):
            self.currentResData = self.resDatas[i]
            startCoeffs = np.array([1.0]*(self.polyOrder+1))
            coeffs.append(self._doFits2(startCoeffs, self.energyRanges[i], True))
        if not self.lockBretWig:
            startCoeffs = coeffs
            coeffs = []
            for i in range(len(self.energyRanges)):
                self.currentResData = self.resDatas[i]
                a = [startCoeffs[i][0].real,startCoeffs[i][0].imag] + startCoeffs[i][1]
                coeffs.append(self._doFits2(a, self.energyRanges[i], False))            
        return coeffs
        
    def _doFits2(self, startCoeffs, eRange, lockBretWig):
        self.lockingBretWig = lockBretWig
        eRange = np.array(eRange)
        epRange = np.array([tw.trace(self.ratEPhaseMats[ene]) for ene in eRange])
        args=(eRange, epRange)
        if self.lockingBretWig:
            coeff = self._leastsq(startCoeffs, args) #Faster method but does not perform well if floating Bret-Wigner parameters.
            return (self.currentResData, coeff.tolist())
        else:
            coeff = self._minimize(startCoeffs, args)
            return (coeff[0]+coeff[1]*1.0j, coeff[2:].tolist())

    def _leastsq(self, startCoeffs, args):
        ret = sp.optimize.leastsq(self._bretWigRes_leastsq, startCoeffs, args)
        print ret
        return ret[0]
    
    def _bretWigRes_leastsq(self, coeffs, ene, ep):
        ret1 = self._bretWigRes(coeffs, ene, ep)
        ret2 = np.sqrt(ret1.real**2 + ret1.imag**2)
        return ret2
    
    def _minimize(self, startCoeffs, args):
        ret = sp.optimize.minimize(self._bretWigRes_minimize, startCoeffs, args, method="Nelder-Mead", options={"maxiter":len(startCoeffs)*2000})
        print ret
        return ret.x
    
    def _bretWigRes_minimize(self, coeffs, ene, ep):
        s = 0.0
        for i in range(len(ene)):
            s += self._bretWigRes(coeffs, ene[i], ep[i]).real**2.0
        return s
    
    def _bretWigRes(self, coeffs, ene, ep):
        if self.lockingBretWig:
            ret = ep - bretWig(self.currentResData, coeffs, ene)
        else:
            ret = ep - bretWig(coeffs[0]+coeffs[1]*1.0j, coeffs[2:], ene)
        return ret
