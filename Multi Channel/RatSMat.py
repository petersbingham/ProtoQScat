# -*- coding: utf-8 -*-
import numpy as np
import scipy.sparse.linalg as sp_sparse_linalg
import sympy as sy
from sympy.matrices import Matrix as sy_matrix
import sympy.polys as sy_polys
import sympy.mpmath as sy_mp
import mpmath
mpmath.mp.dps = 100
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

def checkNumpyVersionForFiles():
  npVer = [int(val) for val in np.__version__.split('.')]
  if npVer[0]>1 or (npVer[0]==1 and npVer[1]>6) or (npVer[0]==1 and npVer[1]==6 and npVer[2]>1): #saveTxt not supported prior.
      return True
  return False


SOLVE_METHOD = 0
class CoeffSolve():
  def __init__(self, numPolyTerms, fitSize, numChannels, verbose):
    if SOLVE_METHOD in [0,1,3,4,5,6,7,8]:
      self.matrixType = 0
    elif SOLVE_METHOD in [2]:
      self.matrixType = 1
    elif SOLVE_METHOD in []:
      self.matrixType = 2
      
    if SOLVE_METHOD in [0,1,8]:
      self.indexType = 0
    elif SOLVE_METHOD in [3,4,5,6,7]:
      self.indexType = 1
    elif SOLVE_METHOD in [2]:
      self.indexType = 2
      
    self.numPolyTerms = numPolyTerms
    self.fitSize = fitSize
    self.numChannels = numChannels 
    self.verbose = verbose
  
  def initialiseMatrices(self):
    if self.matrixType == 0:
      self.sysMat = np.matrix([[0.0]*2*self.numPolyTerms*self.numChannels]*self.fitSize*self.numChannels, dtype=np.complex128)
      self.resVec = np.matrix([[0.0]]*self.fitSize*self.numChannels, dtype=np.complex128)
    elif self.matrixType == 1:
      self.sysMat = mpmath.matrix([[0.0+0.0j]*2*self.numPolyTerms*self.numChannels]*self.fitSize*self.numChannels)       
      self.resVec = mpmath.matrix([[0.0+0.0j]]*self.fitSize*self.numChannels)
    else:
      self.sysMat = sy_matrix([[0.0]*2*self.numPolyTerms*self.numChannels]*self.fitSize*self.numChannels)       
      self.resVec = sy_matrix([[0.0]]*self.fitSize*self.numChannels) 
        
  def setSysElement(self, row, col, val):
    self.sysMat[row,col] = val
  
  def setResult(self, row, col, val):
    self.resVec[row,col] = val 
  
  def solve(self):
    if SOLVE_METHOD == 0:
        self.coeffVec = np.linalg.solve(self.sysMat, self.resVec)
        self._calStr("numpy solve")
    elif SOLVE_METHOD == 1:
        self.coeffVec = np.linalg.lstsq(self.sysMat, self.resVec)[0]
        self._calStr("numpy lstsq")
    elif SOLVE_METHOD == 2:
        self.coeffVec = mpmath.qr_solve(self.sysMat, self.resVec)
        self._calStr("mpmath qr_solve")
    elif SOLVE_METHOD == 3:
        self.coeffVec = self._sparseRet(sp_sparse_linalg.bicg(self.sysMat, self.resVec),"numpy sparse bicg") #, tol=1e-04, maxiter=100000000
    elif SOLVE_METHOD == 4:
        self.coeffVec = self._sparseRet(sp_sparse_linalg.bicgstab(self.sysMat, self.resVec),"numpy sparse bicgstab")
    elif SOLVE_METHOD == 5:
        self.coeffVec = self._sparseRet(sp_sparse_linalg.lgmres(self.sysMat, self.resVec),"numpy sparse lgmres")
    elif SOLVE_METHOD == 6:
        self.coeffVec = self._sparseRet(sp_sparse_linalg.minres(self.sysMat, self.resVec),"numpy sparse minres")
    elif SOLVE_METHOD == 7:
        self.coeffVec = self._sparseRet(sp_sparse_linalg.qmr(self.sysMat, self.resVec),"numpy sparse qmr")
    elif SOLVE_METHOD == 8:
        Q,R = np.linalg.qr(self.sysMat)
        y = np.dot(Q.T,self.resVec)
        self.coeffVec = np.linalg.solve(R,y) 
        self._calStr("numpy qr solve")

  def _calStr(self, typeStr):
    if self.verbose:
        print "Coeffs calculated using " + typeStr 

  def _sparseRet(self, ret, typeStr):
    if self.verbose:
        print "Coeffs calculated using " + typeStr + ". Ret: " + str(ret[1])
    return ret[0]

  def getValue(self, row):
    if self.indexType == 0:
        ret = self.coeffVec[row,0]
    elif self.indexType == 1:
        ret = self.coeffVec[row]
    else:
        ret = self.coeffVec[0][row]
    return ret

  def printSysMatToFile(self):
    if self.matrixType == 0:
        if checkNumpyVersionForFiles():
          np.savetxt("systemMatrix.csv", self.sysMat, delimiter=",")
    else:
        with open("systemMatrix.txt", 'w') as f:
            f.write(str(self.resVec))

class RatSMat(sm.mat):
  def __init__(self, sMatData, kCal, fitSize=None, fitName=None, suppressCmdOut=False, verbose=False):
    self.sMatData = sMatData
    self.suppressCmdOut = suppressCmdOut
    self.verbose = verbose
    self._initData(fitSize)
    self.kCal = kCal
    self.type = TYPE_S
    self.hasCoeffs = False
    self.ene = None
    self._initialiseMatrices()
    self.fitName = fitName
    
    read = False
    if self._doFilesExist() and not ALWAYS_CALCULATE and checkNumpyVersionForFiles():
      try:
        self._readCoefficients()
        read = True
      except Exception as e:
        print "Error reading coefficients will attempt to calculate"

    if not read:
      self._calculateCoefficients()
      if self.fitName and checkNumpyVersionForFiles():
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
    self.coeffSolve = CoeffSolve(self.numPolyTerms, self.fitSize, self.numChannels, self.verbose)
  
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
        self.coeffSolve.initialiseMatrices()
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
                self.coeffSolve.setSysElement(self._row(m,ei), self._alphaIndex(j,ti), alphaCoeff)
                self.coeffSolve.setSysElement(self._row(m,ei), self._betaIndex(j,ti), betaCoeff)
            self.coeffSolve.setResult(self._row(m,ei),0,self._result(m, n, ene))
            ei += 1
        self.coeffSolve.printSysMatToFile()
        self.coeffSolve.solve()
        
        coeffVec = self.coeffSolve.solve()
        self._copyColumnCoeffs(fit, n)
      self._print("Calculated Fit: " + str(fit))
    self.hasCoeffs = True
     
  def _row(self, m, ei):
    return m*self.fitSize + ei
  
  def _alphaIndex(self, m, ti):
    return m*self.numPolyTerms + ti
  
  def _betaIndex(self, m, ti):
    return self.numPolyTerms*self.numChannels + m*self.numPolyTerms + ti
  
  def _primaryAlpha(self, m, n, ene, exp):
    return self.kCal.kl(n,ene,1.0) / self.kCal.kl(m,ene,1.0) * (self.sMatData[ene][m,m]-1.0) * pow(ene,exp)
  
  def _primaryBeta(self, m, n, ene, exp):
    return -1.0j * self.kCal.kl(m,ene,0.0) * self.kCal.kl(n,ene,1.0) * (self.sMatData[ene][m,m]+1.0) * pow(ene,exp)
  
  def _secondaryAlpha(self, m, n, j, ene, exp):
    return self.kCal.kl(n,ene,1.0) / self.kCal.kl(j,ene,1.0) * self.sMatData[ene][m,j] * pow(ene,exp)
  
  def _secondaryBeta(self, m, n, j, ene, exp):
    return -1.0j * self.kCal.kl(j,ene,0.0) * self.kCal.kl(n,ene,1.0) * self.sMatData[ene][m,j] * pow(ene,exp)
    
  def _result(self, m, n, ene):
    num = 0.0
    if m==n:
      num = 1.0
    return num - self.sMatData[ene][m,n]

  def _copyColumnCoeffs(self, fit, n):
    eneKey = self.enes[fit*self.fitSize]
    for ci in range(self.numCoeffs):
      ti = ci-1
      for m in range(self.numChannels):
        if ci==0:
          if m==n:
            self.alphas[eneKey][ci][m,n] = 1.0
        else:
          self.alphas[eneKey][ci][m,n] = complex(self.coeffSolve.getValue(self._alphaIndex(m,ti)))
          self.betas[eneKey][ci][m,n] = complex(self.coeffSolve.getValue(self._betaIndex(m,ti)))
          
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
          t1 = self.kCal.kl(n,self.ene,1.0)/self.kCal.kl(m,self.ene,1.0)*A
          t2 = 1.0j*self.kCal.kl(m,self.ene,0.0)*self.kCal.kl(n,self.ene,1.0)*B
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
  
  def findRoot(self, startingEne, multipler=1.0):
    self.setType(TYPE_FIN)
    try:
        return self._findRootAttempt((startingEne,startingEne+0.01,startingEne+0.02), multipler)
    except ValueError:
        return None
  
  def findRoot_Multi(self, startingEne, multipler=1.0):
    self.setType(TYPE_FIN)
    root = None
    for i in reversed(range(-7,-1)):
        root = self._findRootAttempt((startingEne,startingEne+pow(10,float(i)),startingEne+2.0*pow(10,float(i))), multipler)
        if root is not None:
            break
        root = self._findRootAttempt((startingEne,startingEne-pow(10,float(i)),startingEne-2.0*pow(10,float(i))), multipler)
        if root is not None:
            break
    return root

  def _findRootAttempt(self, startingEne, multipler):
    try:
        return complex(sy_mp.findroot(lambda e: self._getDet(e, multipler), startingEne, solver='muller', maxsteps=1000, tol=2.16840434497100886801e-19))
    except ValueError:
        return None
          
  def findConjRoots(self, startingEne, multipler=1.0):
    return [self.findRoot(startingEne, multipler), self._findRoot(startingEne.real - 1.0j*startingEne.imag, multipler)]

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
            k = sy.symbols('k')
            matLst = []
            for m in range(self.numChannels):
                matLst.append([])
                for n in range(self.numChannels):
                    lm = self.kCal.l(m)
                    ln = self.kCal.l(n)
                    val = 0.0
                    for ci in range(self.numCoeffs):
                        A = alphas[ci][m,n]
                        B = betas[ci][m,n]
                        val += (1.0/2.0)*(1.0/kConversionFactor)**(ci) * ( A*k**(ln-lm+2*ci) - 1.0j*B*k**(ln+lm+1+2*ci) )
                    matLst[len(matLst)-1].append(val)
            deter = sy_matrix(matLst).det()
            roots = self._getRoots_numpy(deter, k)
            if convertToEne:
                mappedRoots = map(lambda val: complex((1.0/kConversionFactor)*val**2), roots)
            else:
                mappedRoots = map(lambda val: complex(val), roots)
            allRoots.extend(mappedRoots)
    return allRoots

  def _getRoots_numpy(self, deter, k):
    coeffs = sy_polys.Poly(deter, k).all_coeffs()
    mappedCoeffs = map(lambda val: complex(val), coeffs)
    return np.roots(mappedCoeffs)     

  def _getRoots_sympy(self, deter, k):
    return sy_polys.Poly(deter, k).nroots(n=25, maxsteps=500, cleanup=True)
