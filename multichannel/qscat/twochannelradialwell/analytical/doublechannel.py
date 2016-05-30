import math
import numpy as np
import scipy.linalg as la
import sympy.mpmath as mpm

import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../../Utilities')
import scattering.matrices as sm
import general.numerical as num
import general.simpplot3D as plot3D
from general.qstype import *

EQUIVALENT_TESTS = False
LIN_ALGEBRA = False
ENEFACTOR = sm.EFROMK_HARTREES

gu = num.Compare()

class DCException(Exception):
  def __init__(self, string):
    self.string = string
  def __str__(self):
    return "DC Error: " + self.string
   
    
class Mats:
  def __init__(self, v1, v2, lam, kCal):
    self.K = Mats.Kmat(kCal)
    self.V = Mats.Vmat(v1, v2, lam, kCal.eneFactor)
    self.A = Mats.Amat(self.K, self.V)
    if LIN_ALGEBRA:
        self.aSq = Mats.aSqMat(self.A)
        self.a = Mats.aMat(self.aSq)
  
  def setEnergy(self, ene):
      self.K.setEnergy(ene)
      if LIN_ALGEBRA:
          self.aSq.calculate()
          self.a.calculate()
  
  def printMats(self):
      print "\nK:"
      print str(self.K) + "\n\nV:"
      print str(self.V) + "\n\nA:"
      print str(self.A) + "\n\nsqrt(A):"
      print str(self.aSq) + "\n\na:"
      print str(self.a)
    
  class Kmat(sm.mat):
    def __init__(self, kCal):
      sm.mat.__init__(self, 2, num.PRECISION)
      self.ene = 0
      self.kCal = kCal
    def setEnergy(self, ene):
      self.ene = ene
    def _getRow(self, i):
      if i==0:
        return [self.k(0), 0]
      else:
        return [0, self.k(1)]
    def k(self, ch):
      return self.kCal.k(ch, self.ene)
    def getPhase(self, ch):
      return self.kCal.getPhase(ch, self.ene)
  
  class Vmat(sm.mat):
    def __init__(self, v1, v2, lam, eneFactor):
      sm.mat.__init__(self, 2, num.PRECISION)
      self.v1 = v1
      self.v2 = v2
      self.lam = lam
      self.eneFactor = eneFactor
    def _getRow(self, i):
      if i==0:
        return [-self.v1*self.eneFactor, -0.5*self.lam*self.eneFactor]
      else:
        return [-0.5*self.lam*self.eneFactor, -self.v2*self.eneFactor]
        
  class Amat(sm.mat):
    def __init__(self, K, V):
      sm.mat.__init__(self, 2, num.PRECISION)
      self.K = K
      self.V = V
    def _getRow(self, i):
      return [QSpow(self.K[i][0],2)-self.V[i][0], QSpow(self.K[i][1],2)-self.V[i][1]]
      
  class aSqMat(sm.mat):
    def __init__(self, A):
      sm.mat.__init__(self, 2, num.PRECISION)
      self.A = A
      self.calculate()
    def calculate(self):
      m = self.A.getMatrix()
      eigvals,eigvecs = np.linalg.eig(m)
      self.aSq = np.diagflat(eigvals) 
    def _getRow(self, i):
      return [self.aSq[i,0], self.aSq[i,1]]
      
  class aMat(sm.mat):
    def __init__(self, aSq):
      sm.mat.__init__(self, 2, num.PRECISION)
      self.aSq = aSq
      self.calculate()
    def calculate(self):
      m = self.aSq.getMatrix()
      self.a = la.sqrtm(m)
    def _getRow(self, i):
      return [self.a[i,0], self.a[i,1]]

class Smat(sm.mat):
  def __init__(self, r0, mats):
    sm.mat.__init__(self, 2, num.PRECISION)
    self.numChannels = 2
    self.mats = mats
    self.r0 = r0

  def setEnergy(self, ene, calUniOp=True):
    self.mats.setEnergy(ene)
    if calUniOp:
      self.uniOpSmat = Smat.UniOpSmat(self)
    else:
      self.uniOpSmat = None
      
  def _getRow(self, i):
    if i==0:
      return [self._S_11(), self._S_12()]
    else:
      return [self._S_21(), self._S_22()]
  
  def findRoot(self, start):
    def eneEqu(self, ene):
      self.setEnergy(ene, False)
      return self._denum(False)
    return mpm.findroot(lambda ene: eneEqu(self, ene), start)
  
  def plotPoles(self, Rs, Is, real):
    self.i = 0
    @np.vectorize
    def eneEqu(self, R, I):
      print self.i
      self.i += 1
      self.setEnergy(R+I*1.0j, False)
      return 1.0 / self._denum(False)
    plot3D.plot(Rs, Is, lambda R, I: eneEqu(self, R, I), real, "Energy", "1/Denum")
      
  #######
      
  def abs(self):
    return Smat.AbsSmat(self)
     
  def uniOp(self):
    return self.uniOpSmat
    
  def isUnitary(self):
    if not gu.complexCompare(self.uniOpSmat[0][0],1.0):
      return False
    if not gu.complexCompare(self.uniOpSmat[0][1],0.0):
      return False
    if not gu.complexCompare(self.uniOpSmat[1][0],0.0):
      return False
    if not gu.complexCompare(self.uniOpSmat[1][1],1.0):
      return False
    return True
  
  class AbsSmat(sm.mat):
    def __init__(self, Smat):
      sm.mat.__init__(self, 2, num.PRECISION)
      self.Smat = Smat
    def _getRow(self, i):
      if i==0:
        return [abs(self.Smat[0][0]), abs(self.Smat[0][1])]
      else:
        return [abs(self.Smat[1][0]), abs(self.Smat[1][1])]
  
  class UniOpSmat(sm.mat):
    def __init__(self, Smat):
      sm.mat.__init__(self, 2, num.PRECISION)
      self.Smat = Smat
      self.calculate()
    def calculate(self):
      m1 = self.Smat.getMatrix()
      m2 = m1.transpose().conjugate()
      self.uniOpS = m1 * m2
    def _getRow(self, i):
      return [self.uniOpS[i,0], self.uniOpS[i,1]]
     
  #######
  
  def _S_11(self):
    return self._g(-self._rho_1(), self._rho_2()) / self._denum() * self._exp(2.0*self._rho_1())
  
  def _S_12(self):
    return 2.0 * (self._zeta_1()-self._zeta_2()) * QSsqrt(self._alp_1()*self._alp_2()*self._rho_1()*self._rho_2()) / self._denum() * self._exp(self._rho_1() + self._rho_2())
  
  def _S_21(self):
    return self._S_12()
  
  def _S_22(self):
    return self._g(self._rho_1(), -self._rho_2()) / self._denum() * self._exp(2.0*self._rho_2())
  
  def _g(self, rho_1, rho_2):
    complex1 = rho_1 * (self._zeta_1()*self._alp_1() - self._zeta_2()*self._alp_2())
    complex2 = rho_2 * (self._zeta_2()*self._alp_1() - self._zeta_1()*self._alp_2())
    real = (self._alp_1()-self._alp_2()) * (rho_1*rho_2 - self._zeta_1()*self._zeta_2())
    return real + (complex1+complex2)*1.0j
  
  def _denum(self, test=True):
    value = self._g(self._rho_1(), self._rho_2())
    if test and gu.complexCompare(value, 0.0):
      raise DCException("_denum: Zero")
    return value
  
  def _exp(self, rho):
    return QSexp(-1.0j*rho)
  
  #######
  
  def _rho_1(self):
    return self._rho_alp(0)
    
  def _rho_2(self):
    return self._rho_alp(1)
  
  def _alp_1(self):
    cal1 = QSpow(self._e_n(0), 2) - QSpow(self._R_alp(0), 2)
    cal2 = QSpow(self._R_alp(1), 2) - QSpow(self._e_n(1), 2)
    if EQUIVALENT_TESTS:
      if not gu.complexCompare(cal1, cal2):
        raise DCException("_alp_1: " + str(cal1) + "   " + str(cal2))
    return cal1
    
  def _alp_2(self):
    cal1 = QSpow(self._e_n(1), 2) - QSpow(self._R_alp(0), 2)
    cal2 = QSpow(self._R_alp(1), 2) - QSpow(self._e_n(0), 2)
    if EQUIVALENT_TESTS:
      if not gu.complexCompare(cal1, cal2):
        raise DCException("_alp_2: " + str(cal1) + "   " + str(cal2))
    return cal1
    
  def _zeta_1(self):
    return self._e_n(0) / QStan(self._e_n(0))
    
  def _zeta_2(self):
    return self._e_n(1) / QStan(self._e_n(1))
    
  #######
  
  def _e_n(self, n):
    if LIN_ALGEBRA:
        cal1 = self.mats.a[n][n] * self.r0
    cal2 = QSsqrt(self._e_n_Sq_alt(n))
    if EQUIVALENT_TESTS:
      if not gu.complexCompare(cal1, cal2):
        raise DCException("_e_n: " + str(cal1) + "   " + str(cal2))
    if LIN_ALGEBRA:
      return cal1
    else:
      return cal2
  
  def _e_n_Sq_alt(self, n):
    first = (QSpow(self._R_alp(0),2.0)+QSpow(self._R_alp(1),2.0)) / 2.0
    second = QSsqrt(  QSpow(QSpow(self._R_alp(0),2.0)-QSpow(self._R_alp(1),2.0),2.0) + 4.0*QSpow(self.mats.V[0][1],2.0)*QSpow(self.r0,4.0)   ) / 2.0
    if n==0:
      return first+second
    else:
      return first-second
  
  def _rho_alp(self, ch):
    return self.mats.K.k(ch) * self.r0
    
  def _R_alp(self, ch):
    cal1 = QSsqrt( self.mats.A[ch][ch] ) * self.r0
    cal2 = QSsqrt( QSpow(self._rho_alp(ch),2) - self.mats.V[ch][ch]*QSpow(self.r0,2.0) )
    if EQUIVALENT_TESTS:
      if not gu.complexCompare(cal1, cal2):
        raise DCException("_R_alp: " + str(cal1) + "   " + str(cal2))
    return cal2
    
  #######



        