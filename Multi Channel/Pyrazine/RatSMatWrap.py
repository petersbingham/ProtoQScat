import numpy as np
import cmath
from Elastic3ChanReader import *
import Scattering.Conversions as conv
import Scattering.Stran as S
from RatSMat import *
from Pyrazine import *

M_NORM = 0
M_PIECEWISE = 1

class PyXSmat(S.XSmat):
    def _getElement(self, m, n):
        self.sourceMat.setEnergy(self.ene)
        t = self.sourceMat[m][n]
        k = self.kCal.k(n,self.ene)
        XS = cmath.pi / pow(k,2.0) * pow(abs(t),2.0)
        return XS * conv.BOHRSQ_to_ANGSQ

def getTotalXS(XSmat):
    XS = 0.0
    for x in np.nditer(XSmat.getMatrix(), flags=['refs_ok']):
        XS += x
    return np.matrix([[XS]])

class RatSMatWrap:
  def __init__(self, fileName, N=None, startIndex=None, endIndex=None, kfitSigns=None, ksigns=None, suppressCmdOut=False):
    if N == -1:
        N = None
    if startIndex == -1:
        startIndex = None 
    if endIndex == -1:
        endIndex = None
    self.numChannels = 3
    self.ene = None
    if kfitSigns is None:
        self.kFitCal = sm.kCalculator(THRESHOLDS, LS, eneFactor=ENEFACTOR)
    else:
        self.kFitCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=kfitSigns, eneFactor=ENEFACTOR)
    if ksigns is None:
        self.kCal = self.kFitCal
    else:
        self.kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=ksigns, eneFactor=ENEFACTOR)
    self.kmats = readkMats(fileName)
    self.mode = M_PIECEWISE
    self.N = N
    self.suppressCmdOut = suppressCmdOut
    if startIndex is not None and endIndex is not None and N is not None:
      self.mode = M_NORM
      self._decimate(startIndex, endIndex, N)
    elif len(self.kmats)%2 != 0:
      self.kmats.pop(ene)
    self.fitName = getFitName(self.kFitCal, startIndex, endIndex)
  
  def _decimate(self, startIndex, endIndex, N):
    self.kmats, step, actualEndIndex, startEne, endEne = sm.decimate(self.kmats, startIndex, endIndex, N)
    if not self.suppressCmdOut:
        print "Decimation:"
        print "  N=%d, Emin=%d(%f), Emax=%d(%f), step=%d" % (N,startIndex,startEne,actualEndIndex,endEne,step)

  def setEnergy(self, ene):
    self.ene = ene
  
  def getMatrix(self):
    if self.ene is None:
      return np.identity(NUMCHANNELS)
    else:
      return self._getSfromKmatrix(self.ene)
  
  def getDiscreteXS(self, title=None, colourCycle=None):
    XSmats = sm.matSequence(title, colourCycle)
    for ene in sorted(self.kmats, key=lambda val: val.real):
        xs = PyXSmat(S.Tmat(self), self.kFitCal)
        xs.setEnergy(ene)
        XSmats[ene] = xs.getMatrix()
    return XSmats
  
  def getTotalDiscreteXS(self, title=None, colourCycle=None):
    totXSmats = sm.matSequence(title, colourCycle)
    xs = PyXSmat(S.Tmat(self), self.kFitCal)
    for ene in sorted(self.kmats, key=lambda val: val.real):
        xs.setEnergy(ene)
        totXSmats[ene] = getTotalXS(xs)
    return totXSmats
  
  def _getSfromKmatrices(self):
    return sm.getSfromKmatrices(self.kmats, NUMCHANNELS)
  
  def _getSfromKmatrix(self, ene):
    return sm.getSfromKmatrix(self.kmats, NUMCHANNELS, ene)
  
  def _getRatSmat(self):
    smats = self._getSfromKmatrices()
    ratSmat = RatSMat(smats, self.kFitCal, fitName=self.fitName, fitSize=self._getRatSmatFitSize(), suppressCmdOut=self.suppressCmdOut)
    ratSmat.kCal = self.kCal
    return ratSmat
  
  def getRatXS(self, title=None, colourCycle=None):
    ratSmat = self._getRatSmat()
    ratXSMats = sm.matSequence(title, colourCycle)
    ratXSmat = PyXSmat(S.Tmat(ratSmat), self.kCal)
    for ene in self.kmats:
      ratXSmat.setEnergy(ene)
      ratXSMats[ene] = ratXSmat.getMatrix()
    return ratXSMats
  
  def getTotalRatXS(self, title=None, colourCycle=None, eneStart=None, eneEnd=None, eneComplexOffset=None, eneSteps=None):
    ratSmat = self._getRatSmat()
    ratTotXSMats = sm.matSequence(title, colourCycle)
    ratXSmat = PyXSmat(S.Tmat(ratSmat), self.kCal)
    if eneStart is None and eneEnd is None and eneSteps is None and eneComplexOffset is None:
        eneRange = self.kmats.keys()
    else:
        d = (float(eneEnd) - float(eneStart)) / float(eneSteps)
        eneRange = [eneStart+d*i for i in range(eneSteps+1)]
    for ene in eneRange:
      ratXSmat.setEnergy(ene)
      ratTotXSMats[ene] = getTotalXS(ratXSmat)
    return ratTotXSMats
    
  def findRoot(self, startingEne, multipler=1.0):
    ratSmat = self._getRatSmat()
    return ratSmat.findRoot(startingEne, multipler)
    
  def findConjRoots(self, startingEne, multipler=1.0):
    ratSmat = self._getRatSmat()
    return ratSmat.findConjRoots(startingEne, multipler)
    
  def findPolyRoots(self):
    ratSmat = self._getRatSmat()
    return ratSmat.findPolyRoots(ENEFACTOR)
    
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