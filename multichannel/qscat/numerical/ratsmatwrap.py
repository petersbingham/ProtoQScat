import numpy as np
import cmath

import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')
sys.path.insert(0,os.getcwd()) #We assume that the specific kreader and description (below) will be here.

from sysdesc import *
from ratsmat import *
from ratsmat.decimator import *
from ratsmat.ephasefitter import *
from matreader import *
import scattering.conversions as conv
from scattering.stran import *

M_NORM = 0
M_PIECEWISE = 1

class T_to_XSmol(T_to_XS):
    def _getElement(self, m, n):
        self.sourceMat.setEnergy(self.ene)
        t = self.sourceMat[m][n]
        k = self.kCal.k(n,self.ene)
        XS = cmath.pi / pow(k,2.0) * pow(abs(t),2.0)
        return XS * conv.BOHRSQ_to_ANGSQ

def getTotalXS(T_to_XS):
    XS = QSsumElements(T_to_XS.getMatrix())
    return QSmatrix([[XS]])

class RatSMatWrap:
    def __init__(self, fileName, N=None, startIndex=None, endIndex=None, kfitSigns=None, ksigns=None, suppressCmdOut=False):
        if N == -1:
            N = None
        if startIndex == -1:
            startIndex = None 
        if endIndex == -1:
            endIndex = None
        self.numChannels = NUMCHANNELS
        self.ene = None
        if kfitSigns is None:
            self.kFitCal = sm.kCalculator(THRESHOLDS, LS, massMult=MASSMULT)
        else:
            self.kFitCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=kfitSigns, massMult=MASSMULT)
        if ksigns is None:
            self.kCal = self.kFitCal
        else:
            self.kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=ksigns, massMult=MASSMULT)
        self.kmats = readkMats(fileName)
        self.mode = M_PIECEWISE
        self.N = N
        self.suppressCmdOut = suppressCmdOut
        self.fileHandler = getFileHandler(self.kFitCal, startIndex if startIndex is not None else 0, endIndex if endIndex is not None else len(self.kmats)-1)
        if startIndex is not None and endIndex is not None and N is not None:
            self.mode = M_NORM
            self._decimate(startIndex, endIndex, N)
        else:
            if len(self.kmats)%2 != 0:
                self.kmats.pop(ene)
  
    def _decimate(self, startIndex, endIndex, N):
        decimator = Decimator(startIndex, endIndex, 0, self.fileHandler)
        self.kmats, decStr = decimator.decimate(self.kmats, N)
        if not self.suppressCmdOut:
            print "Decimation:  " + decStr
    
    def setEnergy(self, ene):
        self.ene = ene
    
    def getMatrix(self):
        if self.ene is None:
            return QSidentity(NUMCHANNELS)
        else:
            return self._getSfromKmatrix(self.ene)
    
    def getDiscreteXS(self, title=None, colourCycle=None):
        XSmats = sm.matSequence(title, colourCycle)
        xs = T_to_XSmol(S_to_T(self), self.kFitCal)
        for ene in sorted(self.kmats, key=lambda val: val.real):
            xs.setEnergy(ene)
            XSmats[ene] = xs.getMatrix()
        return XSmats
    
    def getTotalDiscreteXS(self, title=None, colourCycle=None):
        totXSmats = sm.matSequence(title, colourCycle)
        xs = T_to_XSmol(S_to_T(self), self.kFitCal)
        for ene in sorted(self.kmats, key=lambda val: val.real):
            xs.setEnergy(ene)
            totXSmats[ene] = getTotalXS(xs)
        return totXSmats
    
    def getDiscreteEigenSum(self, title=None, colourCycle=None):
        ratEPhaseMats = sm.matSequence(title, colourCycle)
        ratEPhaseMat = S_to_EPhase(self)
        i = 0 
        for ene in sorted(self.kmats, key=lambda val: val.real):
            ratEPhaseMat.setEnergy(ene)
            ratEPhaseMats[ene] = ratEPhaseMat.getMatrix()
            print i
            i+=1
        return ratEPhaseMats
    
    def getEigenSumFit(self, resDatas, polyOrder, eneStart, eneEnd, eneSteps, title=None, colourCycle=None):
        ratEPhaseMats = self.getDiscreteEigenSum(title, colourCycle)
        ePhaseFitter = EPhaseFitter(ratEPhaseMats, resDatas, polyOrder)
        fitCoeffs = ePhaseFitter.getCoeffs()
        
        ePhaseFitPointsLst = [sm.matSequence(title, colourCycle) for i in range(len(resDatas))]
        
        d = (float(eneEnd) - float(eneStart)) / float(eneSteps)
        eneRange = [eneStart+d*i for i in range(eneSteps+1)]
        for ene in eneRange:
            for i in range(len(resDatas)):
                if inFitRange(resDatas[i], ene):
                    ePhaseFitPointsLst[i][ene] = QSmatrix([[bretWig(fitCoeffs[i][0], fitCoeffs[i][1], ene)]])
                    #ePhaseFitPointsLst[i][ene] = QSmatrix([[resonant(resDatas[i], ene)]]) #resonant part only
        
        return (ratEPhaseMats, ePhaseFitPointsLst)
    
    def _getSfromKmatrices(self):
        return sm.getSfromKmatrices(self.kmats, NUMCHANNELS)
    
    def _getSfromKmatrix(self, ene):
        return sm.getSfromKmatrix(self.kmats, NUMCHANNELS, ene)
    
    def _getRatSmat(self):
        smats = self._getSfromKmatrices()
        ratSmat = RatSMat(smats, self.kFitCal, resultFileHandler=self.fileHandler, fitSize=self._getRatSmatFitSize(), suppressCmdOut=self.suppressCmdOut)
        ratSmat.kCal = self.kCal
        return ratSmat
    
    def getRatXS(self, title=None, colourCycle=None):
        ratSmat = self._getRatSmat()
        ratXSMats = sm.matSequence(title, colourCycle)
        ratXSmat = T_to_XSmol(S_to_T(ratSmat), self.kCal)
        for ene in self.kmats:
            ratXSmat.setEnergy(ene)
            ratXSMats[ene] = ratXSmat.getMatrix()
        return ratXSMats
    
    def getTotalRatXS(self, title=None, colourCycle=None, eneStart=None, eneEnd=None, eneComplexOffset=None, eneSteps=None):
        ratSmat = self._getRatSmat()
        ratTotXSMats = sm.matSequence(title, colourCycle)
        ratXSmat = T_to_XSmol(S_to_T(ratSmat), self.kCal)
        if eneStart is None and eneEnd is None and eneSteps is None and eneComplexOffset is None:
            eneRange = self.kmats.keys()
        else:
            d = (float(eneEnd) - float(eneStart)) / float(eneSteps)
            eneRange = [eneStart+d*i for i in range(eneSteps+1)]
        for ene in eneRange:
            ratXSmat.setEnergy(ene)
            ratTotXSMats[ene] = getTotalXS(ratXSmat)
        return ratTotXSMats
    
    def getRatEigenSum(self, title=None, colourCycle=None, eneStart=None, eneEnd=None, eneComplexOffset=None, eneSteps=None):
        ratSmat = self._getRatSmat()
        ratEPhaseMats = sm.matSequence(title, colourCycle)
        ratEPhaseMat = S_to_EPhase(ratSmat)
        if eneStart is None and eneEnd is None and eneSteps is None and eneComplexOffset is None:
            eneRange = self.kmats.keys()
        else:
            d = (float(eneEnd) - float(eneStart)) / float(eneSteps)
            eneRange = [eneStart+d*i for i in range(eneSteps+1)]
        for ene in eneRange:
            ratEPhaseMat.setEnergy(ene)
            ratEPhaseMats[ene] = ratEPhaseMat.getMatrix()
        return ratEPhaseMats
      
    def findERoot(self, startingEne, multipler=1.0):
        ratSmat = self._getRatSmat()
        return ratSmat.findERoot(startingEne, multipler)
    
    def findERoot_Multi(self, startingEne, multipler=1.0):
        ratSmat = self._getRatSmat()
        return ratSmat.findERoot_Multi(startingEne, multipler)
      
    def findERoot_Conj(self, startingEne, multipler=1.0):
        ratSmat = self._getRatSmat()
        return ratSmat.findERoot_Conj(startingEne, multipler)
      
    def findRoots(self):
        ratSmat = self._getRatSmat()
        return ratSmat.findRoots()
      
    def getFinDetRange(self, startEne, endEne, complexOffset, steps):
        ratSmat = self._getRatSmat()
        return ratSmat.getFinDetRange(startEne, endEne, complexOffset, steps)
    
    def _getRatSmatFitSize(self):
        if self.mode==M_NORM:
            return len(self.kmats) 
        else:
            if self.N is not None:
                return self.N
            else:
                return None