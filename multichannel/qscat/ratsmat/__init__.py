# -*- coding: utf-8 -*-

from auxsolvers import *

TYPE_S = 0
TYPE_FIN = 1
COEFFDIR = os.path.dirname(os.path.realpath(__file__)) + "/CoefficientFiles/"

##########################################################
################### Configuration Here ###################
##########################################################

ALWAYS_CALCULATE = False

PYTYPE_COEFF_SOLVE_METHOD = "numpy_solve"      #"numpy_solve""numpy_lstsq""numpy_sparse_bicg""numpy_sparse_bicgstab""numpy_sparse_lgmres""numpy_sparse_minres""numpy_sparse_qmr""numpy_qr"

EXPANDEDDET_ROOTS_FINDTYPE = "sympy_nroots"    #"numpy_roots""sympy_nroots"
EXPANDEDDET_ROOTS_CLEANWIDTH = 10.0**-3        #Set to None to turn off

SINGLEROOT_FINDTYPE = "muller"                 #"muller""secant"

DISPLAY_PRECISION = 8

######################## Minor ###########################

SYMPY_NROOTS_N = DPS
SYMPY_NROOTS_MAXSTEPS = 5000

##########################################################
##########################################################

class RatSMat(sm.mat):
    def __init__(self, sMatData, kCal, fitSize=None, resultFileHandler=None, suppressCmdOut=False, doCalc=True):
        self.sMatData = sMatData
        self.suppressCmdOut = suppressCmdOut
        self._initData1(fitSize)
        self.kCal = kCal
        self.type = TYPE_S
        self.hasCoeffs = False
        self.ene = None

        self.coeffSolve = CoeffSolve(self.suppressCmdOut, QSMODE, PYTYPE_COEFF_SOLVE_METHOD)
        self.rootSolver = SymDetRoots(self.suppressCmdOut, EXPANDEDDET_ROOTS_FINDTYPE, SYMPY_NROOTS_N, SYMPY_NROOTS_MAXSTEPS)
        self.rootCleaner = RootClean(self.suppressCmdOut, EXPANDEDDET_ROOTS_CLEANWIDTH)

        self.resultFileHandler = resultFileHandler
        if self.resultFileHandler is not None:
            self.resultFileHandler.setFitInfo(self.numFits, self.fitSize)
            self.resultFileHandler.setCoeffRoutine(self.coeffSolve.typeStr)
            self.resultFileHandler.setRootFindRoutine(self.rootSolver.typeStr)
            self.resultFileHandler.setCleanRootParameter(self.rootCleaner.typeStr)
        if ALWAYS_CALCULATE or not canCacheCoefficients(): #We need to set the info for later use before setting reference to None.
            self.resultFileHandler = None
        self.coeffSolve.setResultFileHandler(self.resultFileHandler)
        self.rootSolver.setResultFileHandler(self.resultFileHandler)
        self.rootCleaner.setResultFileHandler(self.resultFileHandler)
        self.fitName = None

        if doCalc:
            self.doCalc()

    # Calculate the coefficients.
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


    # Set the energy values can be accessset for the set energy using the sm.mat interface.
    def setEnergy(self, ene):
        #print "ene: " + str(ene)
        self.ene = ene
        self._calculate()
    # Get the determinant of Fin for the energy set above.
    def getFinDet(self):
        if self.type != TYPE_FIN:
            raise GenUtils.MatException("Wrong type set")
        else:
            return QSdet(self.getMatrix())
    # This returns values of Fin across a specified range.
    def getFinDetRange(self, startEne, endEne, complexOffset, steps):
        xs = np.ndarray((steps,), dtype=float)
        ys = np.ndarray((steps,), dtype=float)
        zs = np.ndarray((steps,), dtype=float)
        self._setType(TYPE_FIN)
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


    # This attempts to find a root using a specified a starting value.
    def findERoot(self, startingEne, multipler=1.0):
        self._setType(TYPE_FIN)
        try:
            return self._findERootAttempt((startingEne,startingEne+0.01,startingEne+0.02), multipler, SINGLEROOT_FINDTYPE)
        except ValueError:
            return None
    # This attempts to find a root using a specified a starting value and another at its conjugate.
    def findERoot_Conj(self, startingEne, multipler=1.0):
        return [self.findERoot(startingEne, multipler), self.findERoot(startingEne.real - 1.0j*startingEne.imag, multipler)]
    # This crude experimental function attempts to find a root by varying the starting across a specified range.
    def findERoot_Multi(self, startingEne, multipler=1.0):
        self._setType(TYPE_FIN)
        root = None
        for i in reversed(range(-7,-1)):
            for j in range(0,10):
                root = self._findERootAttempt((startingEne,startingEne+self._getModifier(i,j),startingEne+2.0*self._getModifier(i,j)), multipler, SINGLEROOT_FINDTYPE)
                if root is not None:
                    break
                root = self._findERootAttempt((startingEne,startingEne-self._getModifier(i,j),startingEne-2.0*self._getModifier(i,j)), multipler, SINGLEROOT_FINDTYPE)
                if root is not None:
                    break
            if root is not None:
                break
        return root


    # This returns all of the roots for the elastic case.
    def findElastickRoots(self, convertToEne=True):
        return self._findRoots(convertToEne, self.rootSolver.getRoots)
    # This cleans the roots for the elastic case.
    def cleanElastickRoots(self, roots):
        #Originally thought that bad roots around start of data set but moved data set away from zero and bad roots did not follow.
        #badRoot = self.kCal.fk(sorted(self.sMatData, key=lambda val: val.real)[0])
        badRoot = 0.0
        return self.rootCleaner.cleanRoots(roots, badRoot)


########################################################################   
###################### Initialisation Functions ########################
########################################################################

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


########################################################################
########################## Coefficient Files ###########################
########################################################################

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


########################################################################
##################### Calculation of Coefficients ######################
########################################################################

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


########################################################################
############################### Numeric ################################
######################################################################## 

#### General Energy Calculation ####

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

    def _setType(self, type):
        if type == TYPE_FIN:
            self.type = TYPE_FIN
        else:
            self.type = TYPE_S
        if self.ene:
            self._calculate()
    
    def _getRow(self, m):    
        if self.hasCoeffs:
            return QSgetRow(self.selMat, m)
        else:
            raise GenUtils.MatException("Calculation Error")

#### Single Root Finders ####

    def _getDet(self, e, multipler=1.0):
        self.setEnergy(e)
        val = multipler*self.getFinDet()
        #print str(e) + "\n" + str(val) + "\n"
        return val
    
    def _getModifier(self, i, j):
        mod1 = QSpow(10,float(i))
        mod2 = QSpow(10,float(i-1))
        return mod1 + j*mod2
    
    def _findERootAttempt(self, startingEne, multipler, type):
        try:
            fun = lambda e: self._getDet(e, multipler)
            if type == 'muller':
                return complex(mpmath.findroot(fun, startingEne, solver='muller'))
            else:
                return complex(mpmath.findroot(fun, startingEne[0], solver='secant'))
        except ValueError:
            return None


########################################################################
############################## Symbolic ################################
######################################################################## 

#### Elastic Wavenumber Calculation ####

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