import os
import sys
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../Utilities')
from General.File import *
import General.QSType as QSType

POLE_STATUS_NEW = "NEW"
POLE_STATUS_LOST = "LOST"
POLE_STATUS_REP = "REP"

class Vals(list):
    def __init__(self, N, S, E):
        super(Vals, self).__init__()
        self.N = N
        self.S = S
        self.E = E
class Roots(Vals):
    pass
class Poles(Vals):
    pass

class Val:
  def __init__(self, k, E):
    self.k = k
    self.E = E
class Root(Val):
  def __init__(self, k, E):
    Val.__init__(self, k, E)
class Pole(Val):
  def __init__(self, k, E, isNew, isLost, convRoots):
    Val.__init__(self, k, E)
    self.isNew = isNew
    self.isLost = isLost
    self.convRoots = convRoots
  def getStatus(self):
    if self.isNew:
      return POLE_STATUS_NEW
    elif self.isLost:
      return POLE_STATUS_LOST
    else:
      return POLE_STATUS_REP
  
class ResultsAnalyser:
    def __init__(self, basePath, rootsDir=None, polesDir=None):
        if "COEFFS-mpmath" in basePath:
            QSType.QSMODE = QSType.MODE_MPMATH
        else:
            QSType.QSMODE = QSType.MODE_NORM
        if rootsDir is not None:
            fileBase = basePath+sep()+rootsDir+sep()+"Roots"+sep()
            allRoots = self._parseFiles(fileBase, Roots, Root)    
        if polesDir is not None:
            fileBase = basePath+sep()+rootsDir+sep()+polesDir+sep()
            allRoots = self._parseFiles(fileBase, Poles, Pole)                
                    
    def _parseFiles(self, fileBase, subContainerClass, typeClass):
        containerClass = []
        for fileName in os.listdir(fileBase):
            if fileName.endswith(".dat"):
                N, S, E = self._getFileParameters(fileName)
                subContainer = subContainerClass(N,S,E)
                containerClass.append(subContainer)
                self._extractValues(fileBase+fileName, containerClass, typeClass)
        return containerClass
            
    def _getFileParameters(self, fileName):
        params = fileName.split("_")
        params[2] = params[2].replace(".dat","")
        N = int(params[0].split("=")[1])
        S = int(params[1].split("=")[1])
        E = int(params[2].split("=")[1])
        return N, S, E
    
    def _extractValues(self, path, container, typeClass):
        path = fixPath(path)
        with open(path, 'r') as f:
            first = True
            for line in f:
                if not first:
                    kstr = line[line.find(']=')+2:line.find('i')]+'j'
                    post = ""
                    if typeClass==Pole:
                        post = " "
                    Estr = line[line.rfind(']=')+2:line.rfind('i'+post)]+'j' 
                    print kstr + "\t" + Estr  
                    if typeClass==Root:
                        type = Root(QSType.QScomplex(kstr), QSType.QScomplex(Estr))
                    else:
                        convRootsStrs = line[line.find("with")+5:].split(" ")
                        convRoots = map(lambda x: map(lambda y: int(y), x[2:-1].replace("]","").split("[")), convRootsStrs)
                        print convRoots
                        type = Pole(QSType.QScomplex(kstr), QSType.QScomplex(Estr), POLE_STATUS_NEW in line, POLE_STATUS_LOST in line, convRoots)
                    container.append(type) 
                first = False      
            
r = ResultsAnalyser("E:\\Peter's Documents\\PhD\\Code\\Git\\ProtoQScat\\Multi Channel\\Results\\Two Channel Radial Well_[1.0, 1.0]_[1.0, 1.0]_1.0_2.0_2.0_0.0_0.0_1.0_1.0_8.0_0.0_1000\\COEFFS-numpy solve {}\\SingleFit","ROOTS-sympy {cleanup True, maxsteps 500, n 25}","Poles_incN_0.05")
roots = Roots(1,2,3)