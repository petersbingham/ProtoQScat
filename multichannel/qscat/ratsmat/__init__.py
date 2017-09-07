# -*- coding: utf-8 -*-
import numpy as np
import scipy.sparse.linalg as sp_sparse_linalg
import sympy as sy
from sympy.matrices import Matrix as sy_matrix
import sympy.polys as sy_polys
from sympy import poly
import mpmath
import collections

import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../utilities')
sys.path.insert(0,base+'.')
import scattering.matrices as sm
from general import *
from general.qstype import *
import general.numerical as num

TYPE_S = 0
TYPE_FIN = 1
COEFFDIR = os.path.dirname(os.path.realpath(__file__)) + "/CoefficientFiles/"

##########################################################
################### Configuration Here ###################

ALWAYS_CALCULATE = False

PYTYPE_COEFF_SOLVE_METHOD = "numpy_solve"      #"numpy_solve""numpy_lstsq""numpy_sparse_bicg""numpy_sparse_bicgstab""numpy_sparse_lgmres""numpy_sparse_minres""numpy_sparse_qmr""numpy_qr"

EXPANDEDDET_ROOTS_FINDTYPE = "sympy_nroots"    #"numpy_roots""sympy_nroots"
EXPANDEDDET_ROOTS_CLEANWIDTH = 10.0**-3        #Set to None to turn off
#EXPANDEDDET_ROOTS_CLEANWIDTH = None

SINGLEROOT_FINDTYPE = "muller"                 #"muller""secant"

DISPLAY_PRECISION = 8

######################## Minor ###########################

SYMPY_NROOTS_MAXSTEPS = 5000
SYMPY_NROOTS_N = DPS

##########################################################
##########################################################


def _canCacheCoefficients():
    if QSMODE == MODE_NORM:
        npVer = [int(val) for val in np.__version__.split('.')]
        if npVer[0]>1 or (npVer[0]==1 and npVer[1]>6) or (npVer[0]==1 and npVer[1]==6 and npVer[2]>1): #saveTxt not supported prior.
            return True
    else:
        return True
    return False

class AuxHelper():
    def __init__(self, suppressCmdOut):
        self.printed = False
        self.suppressCmdOut = suppressCmdOut
    
    def setResultFileHandler(self, resultFileHandler):
        self.resultFileHandler = resultFileHandler
    
    def _printCalStr(self, nounStr, verbStr, wereLoaded):
        if not self.suppressCmdOut and not self.printed:
            addStr = ""
            if wereLoaded:
                addStr = "had been "
            print nounStr + " " + addStr + verbStr + " using " + self.typeStr 
            self.printed = True

    def _startLogAction(self, string):
        if self.resultFileHandler:
            self.resultFileHandler.startLogAction(string)

    def _endLogAction(self, string):
        if self.resultFileHandler:
            self.resultFileHandler.endLogAction(string)

#v2: In _findRoots function matrix elements with poly(...): matLst[len(matLst)-1].append(poly(val))
class ExpandedDetSolve(AuxHelper):
    def __init__(self, suppressCmdOut):
        AuxHelper.__init__(self, suppressCmdOut)
        self.sympy_detArgs = {'method':'berkowitz'} #'bareis''berkowitz''det_LU'
        self.sympy_polyArgs = {} #'lex''grlex'
        self.numpy_rootsArgs = {}
        self.sympy_nrootsArgs = {'n':SYMPY_NROOTS_N, 'maxsteps':SYMPY_NROOTS_MAXSTEPS, 'cleanup':True}
        typeStrStart = "v2_sympy_det"+getArgDesc(sy_matrix.det, self.sympy_detArgs)+",sympy_Poly"+getArgDesc(sy_polys.Poly.__new__, self.sympy_polyArgs)
        if EXPANDEDDET_ROOTS_FINDTYPE == "numpy_roots":
            self.typeStr = typeStrStart+",numpy_roots"+getArgDesc(np.roots, self.numpy_rootsArgs)
        elif EXPANDEDDET_ROOTS_FINDTYPE == "sympy_nroots":
            self.typeStr = typeStrStart+",sympy_nroots"+getArgDesc(sy_polys.polytools.nroots, self.sympy_nrootsArgs)

    def getRoots(self, mat, k, **args):
        self._startLogAction("getRoots")
        if EXPANDEDDET_ROOTS_FINDTYPE == "numpy_roots":
            ret = self._getRoots_numpy_roots(mat, k)
        elif EXPANDEDDET_ROOTS_FINDTYPE == "sympy_nroots":
            ret = self._getRoots_sympy_Poly_nroots(mat, k)
        self._endLogAction("getRoots")
        return ret
      
    def _getRoots_numpy_roots(self, mat, k):
        self._startLogAction("mat.det")
        deter = mat.det(**self.sympy_detArgs)
        self._endLogAction("mat.det")
        
        self._startLogAction("sy_polys.Poly")
        a = sy_polys.Poly(deter, k, **self.sympy_polyArgs)
        self._endLogAction("sy_polys.Poly")
        
        self._startLogAction("all_coeffs")
        coeffs = a.all_coeffs()
        self._endLogAction("all_coeffs")
        
        mappedCoeffs = map(lambda val: complex(val), coeffs)
        
        self._startLogAction("np.roots")
        ret = np.roots(mappedCoeffs, **self.numpy_rootsArgs)
        self._endLogAction("np.roots")
        
        self.printCalStr()
        return ret 
    
    def _getRoots_sympy_Poly_nroots(self, mat, k):
        self._startLogAction("mat.det")
        deter = mat.det(**self.sympy_detArgs)
        self._endLogAction("mat.det")
        
        self._startLogAction("sy_polys.Poly")
        a = sy_polys.Poly(deter, k, **self.sympy_polyArgs) 
        self._endLogAction("sy_polys.Poly") 
        
        self._startLogAction("nroots")
        ret = a.nroots(**self.sympy_nrootsArgs)   
        self._endLogAction("nroots") 
        
        self.printCalStr()
        return ret 
    
    def printCalStr(self, wereLoaded=False):
        self._printCalStr("Roots", "calculated", wereLoaded)

class RootClean(AuxHelper):
    def __init__(self, suppressCmdOut):
        AuxHelper.__init__(self, suppressCmdOut)
        if EXPANDEDDET_ROOTS_CLEANWIDTH is not None:
            self.typeStr = "dk" + str(EXPANDEDDET_ROOTS_CLEANWIDTH) + "," + str(EXPANDEDDET_ROOTS_CLEANWIDTH) + "i"
        else:
            self.typeStr = None

    def cleanRoots(self, roots, badRoot):
        self._startLogAction("Cleaning roots")
        cleanRoots = []
        rejectRoots = []
        if self.typeStr is not None:
            cmp = num.Compare(EXPANDEDDET_ROOTS_CLEANWIDTH) 
            for root in roots:
                if not cmp.complexCompare(abs(root), badRoot):
                    cleanRoots.append(root)
                else:
                    rejectRoots.append(root)
        else:
            cleanRoots = roots
        self.printCalStr()
        self._endLogAction("Cleaning roots") 
        return cleanRoots, rejectRoots
        
    def printCalStr(self, wereLoaded=False):
        self._printCalStr("Roots", "cleaned", wereLoaded)

class CoeffSolve(AuxHelper):
    def __init__(self, suppressCmdOut):
        AuxHelper.__init__(self, suppressCmdOut)
        self.printed = False
        self._action(0)
        self.coeffVec = None
        
    def setValues(self, numPolyTerms, fitSize, numChannels):  
        self.numPolyTerms = numPolyTerms
        self.fitSize = fitSize
        self.numChannels = numChannels 
    
    def initialiseMatrices(self):
        if self.matrixType==0 or self.matrixType==1:
            self.sysMat = QSmatrix(self._getSysMatInit())
            self.resVec = QSmatrix(self._getResVecInit())
        else:
            self.sysMat = sy_matrix(self._getSysMatInit())       
            self.resVec = sy_matrix(self._getResVecInit()) 
    
    def _getSysMatInit(self):
        return [[0.0]*2*self.numPolyTerms*self.numChannels]*self.fitSize*self.numChannels
          
    def _getResVecInit(self):
        return [[0.0]]*self.fitSize*self.numChannels
          
    def setSysElement(self, row, col, val):
        self.sysMat[row,col] = val
    
    def setResult(self, row, col, val):
        self.resVec[row,col] = val 
    
    def solve(self):
        self._startLogAction(self.typeStr)
        self._action()
        self._endLogAction(self.typeStr)
    
    def _action(self, act=1):
        if QSMODE == MODE_MPMATH:
            args = {}
            if act==0:
                self.typeStr = "mpmath_qr_solve_dps"+getArgDesc(mpmath.qr_solve, args)+" DPS"+str(DPS)
                self.matrixType = 1
                self.indexType = 1
            else:
                self.coeffVec = mpmath.qr_solve(self.sysMat, self.resVec, **args)
                self.printCalStr()
        else:
            if PYTYPE_COEFF_SOLVE_METHOD == "numpy_solve":
                args = {}
                if act==0:
                    self.typeStr = "numpy_solve"+getArgDesc(np.linalg.solve, args)
                    self.matrixType = 0
                    self.indexType = 0
                else:
                    self.coeffVec = np.linalg.solve(self.sysMat, self.resVec, **args)
                    self.printCalStr()
            elif PYTYPE_COEFF_SOLVE_METHOD == "numpy_lstsq":
                args = {}
                if act==0:
                    self.typeStr = "numpy_lstsq"+getArgDesc(np.linalg.lstsq, args)
                    self.matrixType = 0
                    self.indexType = 0
                else:
                    self.coeffVec = np.linalg.lstsq(self.sysMat, self.resVec, **args)[0]
                    self.printCalStr()
            elif PYTYPE_COEFF_SOLVE_METHOD == "numpy_sparse_bicg": 
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_bicg"+getArgDesc(sp_sparse_linalg.bicg, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.bicg(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=10*len(self.resVec)
            elif PYTYPE_COEFF_SOLVE_METHOD == "numpy_sparse_bicgstab":
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_bicgstab"+getArgDesc(sp_sparse_linalg.bicgstab, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.bicgstab(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=10*len(self.resVec)
            elif PYTYPE_COEFF_SOLVE_METHOD == "numpy_sparse_lgmres":
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_lgmres"+getArgDesc(sp_sparse_linalg.lgmres, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.lgmres(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=1000
            elif PYTYPE_COEFF_SOLVE_METHOD == "numpy_sparse_minres":
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_minres"+getArgDesc(sp_sparse_linalg.minres, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.minres(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=5*self.sysMat.shape[0]
            elif PYTYPE_COEFF_SOLVE_METHOD == "numpy_sparse_qmr":
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_qmr"+getArgDesc(sp_sparse_linalg.qmr, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.qmr(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=10*len(self.resVec)
            elif PYTYPE_COEFF_SOLVE_METHOD == "numpy_qr":
                args_qr = {}
                args_s = {}
                if act==0:
                    self.typeStr = "numpy_qr"+getArgDesc(np.linalg.qr, args_qr)+",numpy_solve"+getArgDesc(np.linalg.solve, args_s)
                    self.matrixType = 0
                    self.indexType = 0
                else:
                    Q,R = np.linalg.qr(self.sysMat, **args_qr)
                    y = np.dot(Q.T,self.resVec)
                    self.coeffVec = np.linalg.solve(R,y, **args_s) 
                    self.printCalStr()
    
    def _sparseRet(self, ret):
        if not self.suppressCmdOut and not self.printed:
            print "Coeffs calculated using " + self.typeStr + ". Ret: " + str(ret[1])
            self.printed = True
        return ret[0]
    
    def getValue(self, row):
        if QSMODE == MODE_MPMATH:
            ret = self.coeffVec[0][row]
        else:
            if self.indexType == 0:
                ret = self.coeffVec[row,0]
            elif self.indexType == 1:
                ret = self.coeffVec[row]
            else:
                ret = self.coeffVec[0][row]
        return ret
    
    def printMatricesToFile(self, fit, n):
        if self.coeffVec is not None and self.resultFileHandler is not None:
            path = self.resultFileHandler.getCoeffFilePath()
            if not os.path.isdir(path):
                os.makedirs(path)
            initStr = path+str(fit)+"_"+str(n)+"_"
            if self.matrixType == 0:
                if _canCacheCoefficients():
                    np.savetxt(initStr+"sysMat.txt", self.sysMat, delimiter=",")
                    np.savetxt(initStr+"resVec.txt", self.resVec, delimiter=",")
                    np.savetxt(initStr+"coeffVec.txt", self.coeffVec, delimiter=",")
            else:
                with open(initStr+"sysMat.txt", 'w') as f:
                    f.write(str(self.sysMat))
                with open(initStr+"resVec.txt", 'w') as f:
                    f.write(str(self.resVec))
                with open(initStr+"coeffVec.txt", 'w') as f:
                    f.write(str(self.coeffVec))
    
    def printCalStr(self, wereLoaded=False):
        self._printCalStr("Coeffs", "calculated", wereLoaded)


class RatSMat(sm.mat):
    def __init__(self, sMatData, kCal, fitSize=None, resultFileHandler=None, suppressCmdOut=False, doCalc=True):
        self.sMatData = sMatData
        self.suppressCmdOut = suppressCmdOut
        self._initData1(fitSize)
        self.kCal = kCal
        self.type = TYPE_S
        self.hasCoeffs = False
        self.ene = None
        
        self.coeffSolve = CoeffSolve(self.suppressCmdOut)
        self.rootSolver = ExpandedDetSolve(self.suppressCmdOut)
        self.rootCleaner = RootClean(self.suppressCmdOut)
        
        self.resultFileHandler = resultFileHandler
        if self.resultFileHandler is not None:
            self.resultFileHandler.setFitInfo(self.numFits, self.fitSize)
            self.resultFileHandler.setCoeffRoutine(self.coeffSolve.typeStr)
            self.resultFileHandler.setRootFindRoutine(self.rootSolver.typeStr)
            self.resultFileHandler.setCleanRootParameter(self.rootCleaner.typeStr)
        if ALWAYS_CALCULATE or not _canCacheCoefficients(): #We need to set the info for later use before setting reference to None.
            self.resultFileHandler = None
        self.coeffSolve.setResultFileHandler(self.resultFileHandler)
        self.rootSolver.setResultFileHandler(self.resultFileHandler)
        self.rootCleaner.setResultFileHandler(self.resultFileHandler)
        self.fitName = None
        
        if doCalc:
            self.doCalc()
      
    def doCalc(self):
        self._initData2()
        self._initialiseMatrices()
        
        read = False
        if self.resultFileHandler and self.resultFileHandler.doCoeffFilesExist():
            try:
                self.coeffSolve.printCalStr(True)
                self._readCoefficients()
                read = True
            except Exception as e:
                print "Error reading coefficients will attempt to calculate"
        
        if not read:
            self._calculateCoefficients()
            if self.resultFileHandler:
                self._writeCoefficients() 
        
        self.lastPrintedEne = None
        sm.mat.__init__(self, self.numChannels, DISPLAY_PRECISION)
    
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
    
    def _initData1(self, fitSize):
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
      
    def _initData2(self):
        matShape = QSshape(self.sMatData.itervalues().next())    
        if matShape[0]>=1 and matShape[0]==matShape[1]:
            self.numChannels = matShape[0]
            for ene in self.sMatData:
                matShape = QSshape(self.sMatData[ene])
                if matShape[0]!=self.numChannels or matShape[0]!=matShape[1]:
                    raise sm.MatException("Bad Input: Inconsistent matrices")
            self.coeffSolve.setValues(self.numPolyTerms, self.fitSize, self.numChannels)
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
            sz = self.enes[fit*self.fitSize]
            self.alphas[sz] = self._initialiseCoefficients()
            self.betas[sz] = self._initialiseCoefficients()
        self.selMat = QSmatrix(self._getZeroListMats())
    
    def _initialiseCoefficients(self):
        coeffs = []
        for i in range(0, self.numCoeffs):
            mat = QSmatrix(self._getZeroListMats())
            coeffs.append(mat)
        return coeffs
    
    def _getZeroListMats(self):
        return [[0.0+0.0j]*self.numChannels]*self.numChannels
    
    ############# Coefficient Files ###############
      
    def _readCoefficients(self):
        for fit in range(self.numFits):
            for ci in range(0, self.numCoeffs):
                self.alphas[self.enes[fit*self.fitSize]][ci] = self._readCoefficientsForFit(fit, ci, "A")
                self.betas[self.enes[fit*self.fitSize]][ci] = self._readCoefficientsForFit(fit, ci, "B")
            self._print("Loaded Fit: " + str(fit))
        self.hasCoeffs = True
    
    def _readCoefficientsForFit(self, fit, ci, typeString):
        fileName = self.resultFileHandler.getCoeffFileName(fit, ci, typeString)
        if QSMODE == MODE_NORM:
            return np.asmatrix(np.loadtxt(fileName, dtype=np.complex128, delimiter=","))
        else:
            f = open( fileName, "r" )
            s1 = f.read()
            l1 = s1.split("\n")
            l2 = [self._splitmpRows(s) for s in l1]
            l3 = [map(lambda s:s.replace("[","").replace("]","").replace("[","").replace(")",""),l) for l in l2]
            return mpmath.matrix(l3)
    
    def _splitmpRows(self, s):
        if "(" in s:
            return s.split("(")[1:]
        else:
            return s.split("  ")
    
    def _writeCoefficients(self):
        for fit in range(self.numFits):
            for ci in range(self.numCoeffs):
                self._writeCoefficientsForFit(fit, ci, "A", self.alphas)
                self._writeCoefficientsForFit(fit, ci, "B", self.betas)
    
    def _writeCoefficientsForFit(self, fit, ci, typeString, matRef):
        fileName = self.resultFileHandler.getCoeffFileName(fit, ci, typeString)
        mat = matRef[self.enes[fit*self.fitSize]][ci]
        if QSMODE == MODE_NORM:
            np.savetxt(fileName, mat, delimiter=",")
            self._fixFile(fileName)
        else:
            with open(fileName, 'w') as f:
                f.write(str(mat))
          
    def _fixFile(self, fileName):
        f1 = open(fileName, 'r')
        f2 = open(fileName + "_temp", 'w')
        for line in f1:
            f2.write(line.replace("+-", '-'))
        f1.close()
        f2.close()
        os.remove(fileName)
        os.rename(fileName + "_temp", fileName)
          
    ############# Calculation of Coefficients ###############
    
    def _calculateCoefficients(self):
        if self.resultFileHandler:
            self.resultFileHandler.startLogAction("_calculateCoefficients")
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
                self.coeffSolve.printMatricesToFile(fit, n)
                self.coeffSolve.solve()
                self._copyColumnCoeffs(fit, n)
            self._print("Calculated Fit: " + str(fit))
        if self.resultFileHandler:
            self.resultFileHandler.endLogAction("_calculateCoefficients")
        self.hasCoeffs = True
       
    def _row(self, m, ei):
        return m*self.fitSize + ei
    
    def _alphaIndex(self, m, ti):
        return m*self.numPolyTerms + ti
    
    def _betaIndex(self, m, ti):
        return self.numPolyTerms*self.numChannels + m*self.numPolyTerms + ti
    
    def _primaryAlpha(self, m, n, ene, exp):
        return self.kCal.kl(n,ene,1.0) / self.kCal.kl(m,ene,1.0) * (self.sMatData[ene][m,m]-1.0) * QSpow(ene,exp)
    
    def _primaryBeta(self, m, n, ene, exp):
        return -1.0j * self.kCal.kl(m,ene,0.0) * self.kCal.kl(n,ene,1.0) * (self.sMatData[ene][m,m]+1.0) * QSpow(ene,exp)
    
    def _secondaryAlpha(self, m, n, j, ene, exp):
        return self.kCal.kl(n,ene,1.0) / self.kCal.kl(j,ene,1.0) * self.sMatData[ene][m,j] * QSpow(ene,exp)
    
    def _secondaryBeta(self, m, n, j, ene, exp):
        return -1.0j * self.kCal.kl(j,ene,0.0) * self.kCal.kl(n,ene,1.0) * self.sMatData[ene][m,j] * QSpow(ene,exp)
      
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
                    self.alphas[eneKey][ci][m,n] = QScomplex(self.coeffSolve.getValue(self._alphaIndex(m,ti)))
                    self.betas[eneKey][ci][m,n] = QScomplex(self.coeffSolve.getValue(self._betaIndex(m,ti)))
            
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
            return QSdet(self.getMatrix())  
    
    def _getRow(self, m):    
        if self.hasCoeffs:
            return QSgetRow(self.selMat, m)
        else:
            raise GenUtils.MatException("Calculation Error")
        
    def _calculate(self):
        if self.hasCoeffs:
            Fin = QSmatrix(self._getZeroListMats())
            Fout = QSmatrix(self._getZeroListMats())
            alphas = self._getAlphaSet()
            betas = self._getBetaSet()
            for m in range(self.numChannels):
                for n in range(self.numChannels):
                    A = 0.0
                    B = 0.0
                    for ci in range(self.numCoeffs):
                        exp = ci
                        A += alphas[ci][m,n] * QSpow(self.ene, exp)
                        B += betas[ci][m,n] * QSpow(self.ene, exp)
                    t1 = self.kCal.kl(n,self.ene,1.0)/self.kCal.kl(m,self.ene,1.0)*A
                    t2 = 1.0j*self.kCal.kl(m,self.ene,0.0)*self.kCal.kl(n,self.ene,1.0)*B
                    Fin[m,n] = (t1-t2) / 2.0
                    Fout[m,n] = (t1+t2) / 2.0
            #print str(self.ene) + " , " + str(Fin[0,0]) + " , " + str(Fin[0,1]) + " , " + str(Fin[1,0]) + " , " + str(Fin[1,1])
            if self.type == TYPE_FIN:
                self.selMat = Fin
            else:
                self.selMat = Fout * QSinvert(Fin)  #S-matrix
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
            return self._findRootAttempt((startingEne,startingEne+0.01,startingEne+0.02), multipler, SINGLEROOT_FINDTYPE)
        except ValueError:
            return None
    
    def findRoot_Multi(self, startingEne, multipler=1.0):
        self.setType(TYPE_FIN)
        root = None
        for i in reversed(range(-7,-1)):
            for j in range(0,10):
                root = self._findRootAttempt((startingEne,startingEne+self._getModifier(i,j),startingEne+2.0*self._getModifier(i,j)), multipler, SINGLEROOT_FINDTYPE)
                if root is not None:
                    break
                root = self._findRootAttempt((startingEne,startingEne-self._getModifier(i,j),startingEne-2.0*self._getModifier(i,j)), multipler, SINGLEROOT_FINDTYPE)
                if root is not None:
                    break
            if root is not None:
                break
        return root
    
    def _getModifier(self, i, j):
        mod1 = QSpow(10,float(i))
        mod2 = QSpow(10,float(i-1))
        return mod1 + j*mod2
    
    def _findRootAttempt(self, startingEne, multipler, type):
        try:
            fun = lambda e: self._getDet(e, multipler)
            if type == 'muller':
                return complex(mpmath.findroot(fun, startingEne, solver='muller', maxsteps=10000, tol=2.16840434497100886801e-19))
            else:
                return complex(mpmath.findroot(fun, startingEne[0], solver='secant'))
        except ValueError:
            return None
            
    def findConjRoots(self, startingEne, multipler=1.0):
        return [self.findRoot(startingEne, multipler), self.findRoot(startingEne.real - 1.0j*startingEne.imag, multipler)]
    
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
    
    def findPolyRoots(self, convertToEne=True):
        return self._findRoots(convertToEne, self.rootSolver.getRoots)

    def _findRoots(self, convertToEne, fun, **args):
        if self.resultFileHandler:
            self.resultFileHandler.startLogAction("_findRoots")
        if self.hasCoeffs:
            allRoots = []
            for eKey in self.alphas:
                alphas = self.alphas[eKey]
                betas = self.betas[eKey]
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
                            val += (1.0/2.0)*(1.0/self.kCal.massMult)**(ci) * (QSToSympy(A)*k**(ln-lm+2*ci) - sy.I*QSToSympy(B)*k**(ln+lm+1+2*ci) )
                        matLst[len(matLst)-1].append(poly(val))
                        #matLst[len(matLst)-1].append(val) #v1
                mat = sy_matrix(matLst)
                roots = fun(mat, k, **args)
                if convertToEne:
                    mappedRoots = map(lambda val: self.kCal.e(val,True), roots)
                else:
                    mappedRoots = map(lambda val: complex(val), roots)
                allRoots.extend(mappedRoots)
        if self.resultFileHandler:
            self.resultFileHandler.startLogAction("_findRoots")
        return allRoots

    def cleanPolyRoots(self, roots):
        badRoot = self.kCal.fk(sorted(self.sMatData, key=lambda val: val.real)[0])
        return self.rootCleaner.cleanRoots(roots, badRoot)