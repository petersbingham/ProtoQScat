import cmath
import matrices as sm
import general.numerical as num
from general.qstype import *

class TranMat(sm.mat):
    def __init__(self, sourceMat): 
        self.sourceMat = sourceMat
        self.tranMat = None
        matShape = self.getShape()   
        if matShape[0] == matShape[1]:
            size = matShape[0]
        else:
            raise sm.MatException("Bad Input: Matrix not square")
        sm.mat.__init__(self, size, num.PRECISION)
    
    def getShape(self):
        if isinstance(self.sourceMat, TranMat):
            return self.sourceMat.getShape()
        else:
            return QSshape(self.sourceMat.getMatrix())
    
    def getMatrix(self):
        return self.tranMat
      
    def _getRow(self, m):
        return QSgetRow(self.tranMat, m)

class S_to_T(TranMat):
    def setEnergy(self, ene):
        self.sourceMat.setEnergy(ene)
        self._calculate()
      
    def _calculate(self):
        self.tranMat = self.sourceMat.getMatrix() - QSidentity(self.sourceMat.numChannels)
    
class S_to_UniOp(TranMat):
    def setEnergy(self, ene):
        self.sourceMat.setEnergy(ene)
        self._calculate()
      
    def _calculate(self):
        m1 = self.sourceMat.getMatrix()
        m2 = m1.transpose().conjugate()
        self.tranMat = m1 * m2
    
class S_to_EPhase(TranMat):
    def setEnergy(self, ene):
        self.sourceMat.setEnergy(ene)
        self._calculate()
      
    def _calculate(self):
        num = QSidentity(self.sourceMat.numChannels) - self.sourceMat.getMatrix()
        denum = QSidentity(self.sourceMat.numChannels) + self.sourceMat.getMatrix()
        K = 1.0j*QSdot(num, QSinvert(denum))
        self.tranMat = QSatanElements(QSdiagonalise(K))
 
class T_to_XS(TranMat):
    def __init__(self, tMat, kCal):
        TranMat.__init__(self, tMat) 
        self.kCal = kCal
      
    def setEnergy(self, ene):
        self.ene = ene
        self._calculate()
      
    def _calculate(self):
        self.tranMat = QSsqZeros(self.sourceMat.sourceMat.numChannels)
        for m in range(self.size):
            for n in range(self.size):
                self.tranMat[m,n] = self._getElement(m,n)
      
    def _getElement(self, m, n):
        self.sourceMat.setEnergy(self.ene)
        t = self.sourceMat[m][n]
        k = self.kCal.k(n,self.ene)
        l = self.kCal.l(n)
        return cmath.pi / pow(k,2.0) * (2*l+1) * pow(abs(t),2.0)