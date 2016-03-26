import numpy as np
import sympy
import sympy.polys as syp
import sympy.mpmath as mpm
from sympy.matrices import Matrix
import collections

import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../Utilities')
import Scattering.Matrices as sm

PRECISION = 8

TYPE_S = 0
TYPE_FIN = 1

COEFFDIR = os.path.dirname(os.path.realpath(__file__)) + "/CoefficientFiles/"

ALWAYS_CALCULATE = False

class RatSMat(sm.mat):
  def __init__(self, sMatData, kFun, fitSize=None, fitName=None, suppressCmdOut=False):
    self.sMatData = sMatData
    self.suppressCmdOut = suppressCmdOut
    self._initData(fitSize)
    self.kFun = kFun
    self.type = TYPE_S
    self.hasCoeffs = False
    self.ene = None
    self._initialiseMatrices()
    self.fitName = fitName
    if self._doFilesExist() and not ALWAYS_CALCULATE:
      read = True
      self._readCoefficients()
    else:
      read = False
      self._calculateCoefficients()
    if not read and self.fitName:
      self._writeCoefficients() 
    self.lastPrintedEne = None
    sm.mat.__init__(self, self.numChannels, PRECISION)
    self._print("")
  
  def _print(self, msg):
    if not self.suppressCmdOut:
        print msg
  
  def _datPrint(self):
    if not self.suppressCmdOut:
        print "Data:"
        print "  numData: " + str(self.numData)
        print "  numFits: " + str(self.numFits)
        print "Fit:"
        print "  fitSize: " + str(self.fitSize)
        print "  numFitPoints: " + str(self.numPolyTerms)
        print "  numCoeffs: " + str(self.numCoeffs)
  
  def _initData(self, fitSize):
    self.numData = len(self.sMatData)
    if fitSize is None:
      self.fitSize = self.numData
    else:
      self.fitSize = fitSize
    self._checkInput()
    self.enes = [key.real for key in sorted(self.sMatData.keys(), key=lambda val: val.real)]
    self.numFits = self.numData/self.fitSize
    self.numPolyTerms = self.fitSize / 2
    self.numCoeffs = self.numPolyTerms + 1
    self._datPrint()
    matShape = self.sMatData.itervalues().next().shape     
    if matShape[0]>=1 and matShape[0]==matShape[1]:
      self.numChannels = matShape[0]
      for ene in self.sMatData:
        matShape = self.sMatData[ene].shape
        if matShape[0]!=self.numChannels or matShape[0]!=matShape[1]:
          raise sm.MatException("Bad Input: Inconsistent matrices")
    else:
      raise sm.MatException("Bad Input: Matrix not square")
  
  def _checkInput(self):
    if self.numData<2:
        raise sm.MatException("Bad Input: Not enough Data")
    if self.fitSize<2:
        raise sm.MatException("Bad Input: Specified fit size too small")
    if self.fitSize>self.numData:
        raise sm.MatException("Bad Input: Specified fitsize larger than number of data")
    if self.numData%self.fitSize!=0:
        raise sm.MatException("Bad Input: Num of data is not a multiple of the fit size")
    if self.numData%2!=0:
        raise sm.MatException("Bad Input: Number of data not even")
    if self.fitSize%2!=0:
        raise sm.MatException("Bad Input: Fit size not even")
  
  def _initialiseMatrices(self):
    self.alphas = collections.OrderedDict()
    self.betas = collections.OrderedDict()
    for fit in range(self.numFits):
      self.alphas[self.enes[fit*self.fitSize]] = self._initialiseCoefficients()
      self.betas[self.enes[fit*self.fitSize]] = self._initialiseCoefficients()
    self.selMat = np.matrix(self._getZeroListMats(), dtype=np.complex128)
  
  def _initialiseCoefficients(self):
    coeffs = []
    for i in range(0, self.numCoeffs):
        coeffs.append(np.matrix(self._getZeroListMats(), dtype=np.complex128))
    return coeffs

  def _getZeroListMats(self):
    return [[0.0]*self.numChannels]*self.numChannels

############# Coefficient Files ###############
    
  def _doFilesExist(self):
    if self.fitName:
      if os.path.isdir(self._getDirName()):
        for file in self._getAfitNames() + self._getBfitNames():
          if not os.path.isfile(file):
            return False
        return True
    return False
  
  def _getAfitNames(self): 
    return self._getFitNames("A")
    
  def _getBfitNames(self): 
    return self._getFitNames("B")
  
  def _readCoefficients(self):
    for fit in range(self.numFits):
      for i in range(0, self.numCoeffs):
        self.alphas[self.enes[fit*self.fitSize]][i] = self._readCoefficientsForFit(fit, i, "A")
        self.betas[self.enes[fit*self.fitSize]][i] = self._readCoefficientsForFit(fit, i, "B")
      self._print("Loaded Fit: " + str(fit))
    self.hasCoeffs = True

  def _readCoefficientsForFit(self, fit, i, typeString):
    return np.asmatrix(np.loadtxt(self._getFitName(fit, i, typeString), dtype=np.complex128, delimiter=","))
  
  def _writeCoefficients(self):
    for fit in range(self.numFits):
      for i in range(self.numCoeffs):
        if not os.path.isdir(self._getDirName()):
          os.makedirs(self._getDirName())
        self._writeCoefficientsForFit(fit, i, "A", self.alphas)
        self._writeCoefficientsForFit(fit, i, "B", self.betas)

  def _writeCoefficientsForFit(self, fit, i, typeString, matRef):
    fitName = self._getFitName(fit, i, typeString)
    np.savetxt(fitName, matRef[self.enes[fit*self.fitSize]][i], delimiter=",")
    self._fixFile(fitName)
        
  def _fixFile(self, fitName):
    f1 = open(fitName, 'r')
    f2 = open(fitName + "_temp", 'w')
    for line in f1:
        f2.write(line.replace("+-", '-'))
    f1.close()
    f2.close()
    os.remove(fitName)
    os.rename(fitName + "_temp", fitName)
    
  def _getDirName(self):
    return COEFFDIR + str(self.fitName) + "_" + str(len(self.sMatData)) + "_" + str(self.numFits) + "_" + str(self.numPolyTerms)
  
  def _baseCoeffFile(self): 
    return self._getDirName() + "/"
  
  def _getFitNames(self, typeString): 
    fitNames = []
    for fit in range(self.numFits):
      for i in range(0, self.numCoeffs):
        fitNames.append(self._getFitName(fit, i, typeString))
    return fitNames
  
  def _getFitName(self, fit, i, typeString):
    return self._baseCoeffFile() + typeString + "_" + str(fit) + "_" + str(i) + ".dat" 
        
############# Calculation of Coefficients ###############
 
  def _calculateCoefficients(self):
    for fit in range(self.numFits):
      for n in range(self.numChannels):
        resVec = np.matrix([[0.0]]*self.fitSize*self.numChannels, dtype=np.complex128)
        sysMat = np.matrix([[0.0]*2*self.numPolyTerms*self.numChannels]*self.fitSize*self.numChannels, dtype=np.complex128)
        #print sysMat.shape
        #print resVec.shape
        #print self.numPolyTerms
        for m in range(self.numChannels):
          ei = 0
          for ene in self.enes[fit*self.fitSize:(fit+1)*self.fitSize]:
            for ti in range(self.numPolyTerms):  #We have two indices ci (coefficient) and ti (term). We know the first term in the poly expansion so self.numCoeffs = self.numPolyTerms + 1 
              exp = ti+1
              for j in range(self.numChannels): 
                if j==m:
                  alphaCoeff = self._primaryAlpha(m, n, ene, exp)
                  betaCoeff = self._primaryBeta(m, n, ene, exp)
                else:
                  alphaCoeff = self._secondaryAlpha(m, n, j, ene, exp)
                  betaCoeff = self._secondaryBeta(m, n, j, ene, exp)
                #print str(self._row(m,ei)) + "\t" + str(self._alphaIndex(j,ti)) + "   " + str(alphaCoeff)
                #print "\t" + str(self._betaIndex(j,ti)) + "   " + str(betaCoeff)
                sysMat[self._row(m,ei), self._alphaIndex(j,ti)] = alphaCoeff
                sysMat[self._row(m,ei), self._betaIndex(j,ti)] = betaCoeff
            resVec[self._row(m,ei),0] = self._result(m, n, ene)
            ei += 1
        self._printToFile(sysMat)
        coeffVec = np.linalg.solve(sysMat, resVec)
        #coeffVec = np.linalg.lstsq(sysMat, resVec)[0]
        #print coeffVec
        self._copyColumnCoeffs(fit, coeffVec, n)
      self._print("Calculated Fit: " + str(fit))
    self.hasCoeffs = True
     
  def _row(self, m, ei):
    return m*self.fitSize + ei
  
  def _alphaIndex(self, m, ti):
    return m*self.numPolyTerms + ti
  
  def _betaIndex(self, m, ti):
    return self.numPolyTerms*self.numChannels + m*self.numPolyTerms + ti
  
  def _primaryAlpha(self, m, n, ene, exp):
    return self.kFun(n,ene,1.0) / self.kFun(m,ene,1.0) * (self.sMatData[ene][m,m]-1.0) * pow(ene,exp)
  
  def _primaryBeta(self, m, n, ene, exp):
    return -1.0j * self.kFun(m,ene,0.0) * self.kFun(n,ene,1.0) * (self.sMatData[ene][m,m]+1.0) * pow(ene,exp)
  
  def _secondaryAlpha(self, m, n, j, ene, exp):
    return self.kFun(n,ene,1.0) / self.kFun(j,ene,1.0) * self.sMatData[ene][m,j] * pow(ene,exp)
  
  def _secondaryBeta(self, m, n, j, ene, exp):
    return -1.0j * self.kFun(j,ene,0.0) * self.kFun(n,ene,1.0) * self.sMatData[ene][m,j] * pow(ene,exp)
    
  def _result(self, m, n, ene):
    num = 0.0
    if m==n:
      num = 1.0
    return num - self.sMatData[ene][m,n]

  def _copyColumnCoeffs(self, fit, coeffVec, n):
    eneKey = self.enes[fit*self.fitSize]
    for ci in range(self.numCoeffs):
      ti = ci-1
      for m in range(self.numChannels):
        if ci==0:
          if m==n:
            self.alphas[eneKey][ci][m,n] = 1.0
        else:
          self.alphas[eneKey][ci][m,n] = coeffVec[self._alphaIndex(m,ti),0]
          self.betas[eneKey][ci][m,n] = coeffVec[self._betaIndex(m,ti),0]

  def _printToFile(self, mat):
    np.savetxt("systemMatrix.csv", mat, delimiter=",")
          
#########################################################  
  
  def setType(self, type):
    if type == TYPE_FIN:
      self.type = TYPE_FIN
    else:
      self.type = TYPE_S
    if self.ene:
      self._calculate()
  
  def setEnergy(self, ene):
    #print "ene: " + str(ene)
    self.ene = ene
    self._calculate()
  
  def getFinDet(self):
    if self.type != TYPE_FIN:
      raise GenUtils.MatException("Wrong type set")
    else:
      ret = np.linalg.det(self.getMatrix())
      #print "det: " + str(ret) + "\n"
      return ret
  
  def _getRow(self, m):    
    if self.hasCoeffs:
      return self.selMat[m].tolist()[0]
    else:
      raise GenUtils.MatException("Calculation Error")
      
  def _calculate(self):
    if self.hasCoeffs:
      Fin = np.matrix(self._getZeroListMats(), dtype=np.complex128)
      Fout = np.matrix(self._getZeroListMats(), dtype=np.complex128)
      alphas = self._getAlphaSet()
      betas = self._getBetaSet()
      for m in range(self.numChannels):
        for n in range(self.numChannels):
          A = 0.0
          B = 0.0
          for ci in range(self.numCoeffs):
            exp = ci
            A += alphas[ci][m,n] * pow(self.ene, exp)
            B += betas[ci][m,n] * pow(self.ene, exp)
          t1 = self.kFun(n,self.ene,1.0)/self.kFun(m,self.ene,1.0)*A
          t2 = 1.0j*self.kFun(m,self.ene,0.0)*self.kFun(n,self.ene,1.0)*B
          Fin[m,n] = (t1-t2) / 2.0
          Fout[m,n] = (t1+t2) / 2.0
      #print str(self.ene) + " , " + str(Fin[0,0]) + " , " + str(Fin[0,1]) + " , " + str(Fin[1,0]) + " , " + str(Fin[1,1])
      if self.type == TYPE_FIN:
        self.selMat = Fin
      else:
        self.selMat = Fout * np.linalg.inv(Fin)  #S-matrix
    else:
      raise GenUtils.MatException("Calculation Error")
    
  def _getAlphaSet(self):
    return self._getCoeffSet(self.alphas)
  
  def _getBetaSet(self):
    return self._getCoeffSet(self.betas)
    
  def _getCoeffSet(self, coeffs):
    ene, coeffSet = self._calCoeffSet(coeffs) 
    #if self.lastPrintedEne is None or self.lastPrintedEne!=ene:
    #    print "Used set: " + str(ene)
    self.lastPrintedEne = ene
    return coeffSet
      
  def _calCoeffSet(self, coeffs):
    lastEne = None
    for ene in coeffs:
      if len(coeffs) == 1:
        return ene, coeffs[ene]                 #only one set
      else:
        if lastEne is not None:
          if ene >= self.ene.real:
            return lastEne, coeffs[lastEne]     #have moved into range
        elif ene >= self.ene.real:
          return ene, coeffs[ene]               #energy is before start point
        lastEne = ene
    return lastEne, coeffs[lastEne]             #energy is in last set or after end point
  
  def _getDet(self, e, multipler=1.0):
    self.setEnergy(e)
    val = multipler*self.getFinDet()
    #print str(e) + "\n" + str(val) + "\n"
    return val
  
  def _findRoot(self, startingEne, multipler):
    self.setType(TYPE_FIN)
    try:
        return complex(mpm.findroot(lambda e: self._getDet(e, multipler), (startingEne,startingEne+.001,startingEne+.002), solver='muller'))
    except ValueError:
        return None
          
  def findRoot(self, startingEne, multipler=1.0):
    return self._findRoot(startingEne, multipler)
          
  def findConjRoots(self, startingEne, multipler=1.0):
    return [self._findRoot(startingEne, multipler), self._findRoot(startingEne.real - 1.0j*startingEne.imag, multipler)]

  def getFinDetRange(self, startEne, endEne, complexOffset, steps):
    xs = np.ndarray((steps,), dtype=float)
    ys = np.ndarray((steps,), dtype=float)
    zs = np.ndarray((steps,), dtype=float)
    self.setType(TYPE_FIN)
    ene = startEne
    dene = (endEne - startEne) / float(steps)
    for i in range(0,steps):
      xs[i] = ene
      self.setEnergy(ene + complexOffset*1.0j)
      finDet = self.getFinDet()
      ys[i] = finDet.real
      zs[i] = finDet.imag
      ene += dene
    return (xs, ys, zs)

  def findPolyRoots(self, kConversionFactor, convertToEne=True): #kConversionFactor for when converting from k to energy. eg 2.0 for Hartrees.
    if self.hasCoeffs:
        allRoots = []
        for eKey in self.alphas:
            alphas = self.alphas[eKey]
            betas = self.betas[eKey]
            Fin = np.matrix(self._getZeroListMats(), dtype=np.complex128)
            k = sympy.symbols('k')
            matLst = []
            for m in range(self.numChannels):
                matLst.append([])
                for n in range(self.numChannels):
                    val = 0
                    for ci in range(self.numCoeffs):
                        A = alphas[ci][m,n]
                        B = betas[ci][m,n]
                        val += (1.0/2.0)*(1.0/kConversionFactor)**(ci) * ( A*k**(2*ci) - 1.0j*B*k**(2*ci+1) )
                    matLst[len(matLst)-1].append(val)
            deter = Matrix(matLst).det()
            coeffs = syp.Poly(deter, k).all_coeffs()
            mappedCoeffs = map(lambda val: complex(val), coeffs)
            roots = np.roots(mappedCoeffs)
            if convertToEne:
                mappedRoots = map(lambda val: complex((1.0/kConversionFactor)*val**2), roots)
            else:
                mappedRoots = map(lambda val: complex(val), roots)
            allRoots.extend(mappedRoots)
    return allRoots