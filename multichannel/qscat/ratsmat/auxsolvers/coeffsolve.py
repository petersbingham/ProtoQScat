from auxhelper import *

class CoeffSolve(AuxHelper):
    def __init__(self, suppressCmdOut, qsmode, solve_method):
        AuxHelper.__init__(self, suppressCmdOut)
        self.printed = False
        self.coeffVec = None
        self.qsmode = qsmode
        self.solve_method = solve_method
        self._action(0)
        
    def setValues(self, numPolyTerms, fitSize, numChannels):  
        self.numPolyTerms = numPolyTerms
        self.fitSize = fitSize
        self.numChannels = numChannels 
    
    def initialiseMatrices(self):
        if self.matrixType==0 or self.matrixType==1:
            self.sysMat = MTmatrix(self._getSysMatInit())
            self.resVec = MTmatrix(self._getResVecInit())
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
        if self.qsmode == MODE_MPMATH:
            args = {}
            if act==0:
                self.typeStr = "mpmath_qr_solve_dps"+getArgDesc(mpmath.qr_solve, args)+", DPS "+str(DPS)
                self.matrixType = 1
                self.indexType = 1
            else:
                self.coeffVec = mpmath.qr_solve(self.sysMat, self.resVec, **args)
                self.printCalStr()
        else:
            if self.solve_method == "numpy_solve":
                args = {}
                if act==0:
                    self.typeStr = "numpy_solve"+getArgDesc(np.linalg.solve, args)
                    self.matrixType = 0
                    self.indexType = 0
                else:
                    self.coeffVec = np.linalg.solve(self.sysMat, self.resVec, **args)
                    self.printCalStr()
            elif self.solve_method == "numpy_lstsq":
                args = {}
                if act==0:
                    self.typeStr = "numpy_lstsq"+getArgDesc(np.linalg.lstsq, args)
                    self.matrixType = 0
                    self.indexType = 0
                else:
                    self.coeffVec = np.linalg.lstsq(self.sysMat, self.resVec, **args)[0]
                    self.printCalStr()
            elif self.solve_method == "numpy_sparse_bicg": 
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_bicg"+getArgDesc(sp_sparse_linalg.bicg, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.bicg(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=10*len(self.resVec)
            elif self.solve_method == "numpy_sparse_bicgstab":
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_bicgstab"+getArgDesc(sp_sparse_linalg.bicgstab, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.bicgstab(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=10*len(self.resVec)
            elif self.solve_method == "numpy_sparse_lgmres":
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_lgmres"+getArgDesc(sp_sparse_linalg.lgmres, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.lgmres(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=1000
            elif self.solve_method == "numpy_sparse_minres":
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_minres"+getArgDesc(sp_sparse_linalg.minres, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.minres(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=5*self.sysMat.shape[0]
            elif self.solve_method == "numpy_sparse_qmr":
                args = {}
                if act==0:
                    self.typeStr = "numpy_sparse_qmr"+getArgDesc(sp_sparse_linalg.qmr, args)
                    self.matrixType = 0
                    self.indexType = 1
                else:
                    self.coeffVec = self._sparseRet(sp_sparse_linalg.qmr(self.sysMat, self.resVec, **args))#, tol=1e-05, maxiter=10*len(self.resVec)
            elif self.solve_method == "numpy_qr":
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
        if self.qsmode == MODE_MPMATH:
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
                if canCacheCoefficients():
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