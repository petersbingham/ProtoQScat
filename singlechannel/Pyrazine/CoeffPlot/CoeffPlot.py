import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')
sys.path.insert(0,base+'/../../../Utilities')

import sys
import csv
import cmath
import numpy as np
import math

import general.simpplot as sp
sp.setImgSize(15,6)
from PyrazineDataVals import *

debugFile1 = open('debug1.txt','w')
debugFile2 = open('debug2.txt','w')

class Calculator:
  def _Xsec(self, k, smat):
    bohrSq = cmath.pi / (pow(k,2)) * pow(abs(1-smat),2.0)
    SI = bohrSq * pow(5.292E-11,2)
    barns = SI / 1E-28
    return bohrSq

  def _phase(self, smat):
      return (cmath.log(smat) / 2.0j).real
    
class KCalculator(Calculator):
  def __init__(self, fileName):
    self.energies = []
    self.ks = []
    with open(fileName, 'rb') as file:
      for i in range(0,5):
        file.readline()
      for line in file:
        nums = line.split()
        if len(nums) == 4:
          self.energies.append(float(nums[3]))
        elif len(nums) == 1:
          self.ks.append(float(nums[0]))
        else:
          raise Exception("Bad line")
    if len(self.energies) != len(self.ks):
      raise Exception("Bad data")
      
  def __len__(self):
    return len(self.energies)
    
  def k(self, index):
    if index >= len(self.energies):
      raise Exception("Outside range")
    return math.sqrt(self.energies[index])
  
  def rydbergs(self, index):
      return self.energies[index]
  
  def eV(self, index):
      return 27.21138505 * self.energies[index]
  
  def Smat(self, index):
    return (1+1j*self.ks[index]) / (1-1j*self.ks[index])
    
  def Xsec(self, index):
    k = self.k(index)
    smat = self.Smat(index)
    return self._Xsec(k, smat)

  def phase(self, index):
    smat = self.Smat(index)
    debugFile1.write(str(smat)+'\n')
    return self._phase(smat)
  
class PadeCalculator(Calculator):
  def __init__(self, fileName, N):
    self.title = None
    self.coeffs = []
    with open(fileName, 'rb') as coeffs:
      foundN = False
      first = True
      for row in coeffs:
        if not foundN:
          try:
            if row[0] == 'N':
              Nstr = row[2:].split(",")[0]
              if int(Nstr) == N:
                self._createTitle(row)
                foundN = True
          except IndexError:
            pass
        else:
          if len(row)>3:
            print complex(row.replace('i','j'))
            if not first:
              self.coeffs.append(complex(row.replace('i','j')))
            else:
              first = False
          else:
            break
  
  def _createTitle(self, row):
    rowSplit = row.split(',')
    self.title = rowSplit[0]
    for val in rowSplit[1:]:
        valSplit = val.split('=')
        self.title += "," + valSplit[0] + "=" + stepToRydStr(valSplit[1])
  
  def rydbergs(self, k):
    return pow(k, 2.0)
  
  def eV(self, k):
    a = pow(k / 5.292E-11 * 1.054571726E-34,2.0)
    b = 9.10938291E-31
    return ( a / b) / 1.602176565E-19
  
  def Smat(self, k):
    num = 1.0
    denum = 1.0
    for n in range(len(self.coeffs)):
      mult = self.coeffs[n] * complex(math.pow(k,n+1), 0)
      num += mult
      denum += math.pow(-1,n+1) * mult
    return num / denum
    
  def Xsec(self, k):
    smat = self.Smat(k)
    return self._Xsec(k, smat)

  def phase(self, k):
    smat = self.Smat(k)
    debugFile2.write(str(smat)+'\n')
    return self._phase(smat)

kc = KCalculator('tempfort.19')
fileName = sys.argv[1]
if len(sys.argv) > 2:
  type = float(sys.argv[2])
else:
  type = 1
if len(sys.argv) > 3:
  N = float(sys.argv[3])
else:
  N = 65
if len(sys.argv) > 4:
  startk = float(sys.argv[4])
else:
  startk = kc.k(0)
if len(sys.argv) > 5:
  endk = float(sys.argv[5])
else:
  endk = kc.k(len(kc)-1)
  
pc = PadeCalculator(fileName, N)
if type==1:
  title = 'Pyrazine(Au) Scattering Cross Sections (single channel) '
  ylabel = 'Cross Section (bohr^2)'
  fileSupp = 'XC'
elif type==2:
  title = 'Pyrazine(Au) Scattering Phase Shifts (single channel) '
  ylabel = 'Phase (rad)'
  fileSupp = 'Phase'
elif type==3:
  title = 'Pyrazine(Au) Real S-matrix (single channel) '
  ylabel = ''
  fileSupp = 'Sreal'
elif type==4:
  title = 'Pyrazine(Au) Imag S-matrix (single channel) '
  ylabel = ''
  fileSupp = 'Simag'

xE1 = np.ndarray((len(kc),), dtype=float)
y1 = np.ndarray((len(kc),), dtype=float)
for i in range(len(kc)):
  xE1[i] = kc.rydbergs(i)
  if type==1:
    y1[i] = kc.Xsec(i)
  elif type==2:
    y1[i] = kc.phase(i)
  elif type==3:
    y1[i] = kc.Smat(i).real
  elif type==4:
    y1[i] = kc.Smat(i).imag
sp.plotSingle2(title, [xE1], [y1], xlabel='Electron Energy (rydbergs)',ylabel=ylabel,legends=["From raw K-matrix data"],path="results/"+fileSupp+"_Raw_"+pc.title+".png")

steps = len(kc)
xk = np.arange(startk, endk, (endk-startk)/steps);
xE2 = np.ndarray((steps,), dtype=float)
y2 = np.ndarray((steps,), dtype=float)
for i in range(0,len(xk)):
  xE2[i] = pc.rydbergs(xk[i])
  if type==1:
    y2[i] = pc.Xsec(xk[i])
  elif type==2:
    y2[i] = pc.phase(xk[i])
  elif type==3:
    y2[i] = pc.Smat(xk[i]).real
  elif type==4:
    y2[i] = pc.Smat(xk[i]).imag
sp.plotSingle2(title+pc.title, [xE2], [y2], xlabel='Electron Energy (rydbergs)',ylabel=ylabel,legends=["From rational S-matrix"],path="results/"+fileSupp+"_Rat_"+pc.title+".png")

sp.plotSingle2(title+pc.title, [xE1,xE2], [y1,y2], xlabel='Electron Energy (rydbergs)',ylabel=ylabel,legends=["From raw K-matrix data","From rational S-matrix"],path="results/"+fileSupp+"_Com_"+pc.title+".png")

debugFile1.close()
debugFile2.close()
