import numpy as np
from Elastic3ChanReader import *
import Scattering.Stran as S
from RatSMat import *

M_NORM = 0
M_PIECEWISE = 1

class RatSMatWrap:
  def __init__(self, fileName, N=None, startIndex=None, endIndex=None, kfitSigns=None, ksigns=None, suppressCmdOut=False):
    self.numChannels = 3
    self.ene = None
    if kfitSigns is None:
        self.kFitCal = sm.kCalculator([0.0,0.0,0.0], eneFactor=ENEFACTOR)
    else:
        self.kFitCal = sm.kCalculator([0.0,0.0,0.0], ktype=sm.K_SIGN, ksigns=kfitSigns, eneFactor=ENEFACTOR)
    if ksigns is None:
        self.kCal = self.kFitCal
    else:
        self.kCal = sm.kCalculator([0.0,0.0,0.0], ktype=sm.K_SIGN, ksigns=ksigns, eneFactor=ENEFACTOR)
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
    self.kmats, step, endIndex = sm.decimate(self.kmats, startIndex, endIndex, N)
    if not self.suppressCmdOut:
        print "Decimation:"
        print "  step: " + str(step)

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
      xs = S.XSmat(S.Tmat(self), self.kFitCal)
      xs.setEnergy(ene)
      XSmats[ene] = xs.getMatrix()
    return XSmats
  
  def _getSfromKmatrices(self):
    return sm.getSfromKmatrices(self.kmats, NUMCHANNELS)
  
  def _getSfromKmatrix(self, ene):
    return sm.getSfromKmatrix(self.kmats, NUMCHANNELS, ene)
  
  def _getRatSmat(self):
    smats = self._getSfromKmatrices()
    ratSmat = RatSMat(smats, self.kFitCal.kl, fitName=self.fitName, fitSize=self._getRatSmatFitSize(), suppressCmdOut=self.suppressCmdOut)
    ratSmat.kFun = self.kCal.kl
    return ratSmat
  
  def getRatXS(self, title=None, colourCycle=None):
    ratSmat = self._getRatSmat()
    ratXSMats = sm.matSequence(title, colourCycle)
    ratXSmat = S.XSmat(S.Tmat(ratSmat), self.kCal)
    for ene in self.kmats:
      ratXSmat.setEnergy(ene)
      ratXSMats[ene] = ratXSmat.getMatrix()
    return ratXSMats
    
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