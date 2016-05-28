import os
import sys
from tabulate import tabulate
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../Utilities')
from General import *
from General.File import *
import General.Numerical as num
import PoleFinder
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
  
DISPLAY_DIFFPRESISION = 16

class ResultsAnalyser:
    def __init__(self, basePath, rootsDir=None, polesDir=None):
        if "COEFFS-mpmath" in basePath:
            QSType.QSMODE = QSType.MODE_MPMATH
        else:
            QSType.QSMODE = QSType.MODE_NORM
        if rootsDir is not None:
            fileBase = basePath+sep()+rootsDir+sep()+"Roots"+sep()
            self.allRoots = self._parseFiles(fileBase, Roots, Root)
        if polesDir is not None:
            fileBase = basePath+sep()+rootsDir+sep()+polesDir+sep()
            self.allPoles = self._parseFiles(fileBase, Poles, Pole)
        self.distFactor = float(polesDir[polesDir.find("_df")+3:polesDir.find("_zero")])
        self.zeroValExp = int(polesDir[polesDir.find("_zero")+6:])
        zeroVal = 10**(self.zeroValExp)
        self.ratCmp = num.RationalCompare(zeroVal, self.distFactor)

    def _parseFiles(self, fileBase, subContainerClass, typeClass):
        containerClass = []
        for fileName in os.listdir(fixPath(fileBase)):
            if fileName.endswith(".dat"):
                N, S, E = self._getFileParameters(fileName)
                subContainer = subContainerClass(N,S,E)
                containerClass.append(subContainer)
                self._extractValues(fileBase+fileName, subContainer, typeClass)
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
                if not first and PoleFinder.COMPLETE_STR not in line:
                    kstr = line[line.find(']=')+2:line.find('i')]+'j'
                    post = ""
                    if typeClass==Pole:
                        post = " "
                    Estr = line[line.rfind(']=')+2:line.rfind('i'+post)]+'j' 
                    #print kstr + "\t" + Estr  
                    if typeClass==Root:
                        type = Root(QSType.QScomplex(kstr), QSType.QScomplex(Estr))
                    else:
                        convRootsStrs = line[line.find("with")+5:].split(" ")
                        convRoots = map(lambda x: map(lambda y: int(y), x[2:-1].replace("]","").split("[")), convRootsStrs)
                        type = Pole(QSType.QScomplex(kstr), QSType.QScomplex(Estr), POLE_STATUS_NEW in line, POLE_STATUS_LOST in line, convRoots)
                    container.append(type) 
                first = False      
                
    def createPoleTable(self):
        poleSets = []
        poleLastValues = []
        for poles in self.allPoles:
            for pole in poles:
                found = False
                for i in range(len(poleLastValues)):
                    poleLastValue = poleLastValues[i]
                    if self.ratCmp.isClose(poleLastValue, pole.k):
                        poleSets[i][poles.N] = pole
                        poleLastValues[i] = pole.k
                        found = True
                        break
                if not found:
                    poleSets.append({poles.N:pole})
                    poleLastValues.append(pole.k)
        table = []
        for poleSet in poleSets:
            table.append([" ", "_", "_", "_"])
            first = True
            for N in sorted(poleSet.keys()):
                pole = poleSet[N]
                if first:
                    self._writePriorRoots(table, pole)
                first = False
                table.append([N, pole.getStatus(), self._v(pole.E.real), self._v(pole.E.imag)])
        outStr = getFormattedHTMLTable(table,'.'+str(self.zeroValExp)+'f', headers=["N","Status", "pole.E.real", "pole.E.imag"])
        with open("PoleAnalysers/out.txt", 'w') as f:
            f.write(outStr)
        print outStr
    
    def _writePriorRoots(self, table, initPole):
        for priorRoot in initPole.convRoots[1:]:
            N = priorRoot[0]
            I = priorRoot[1]
            Nroots = self._getRootsForN(N)
            table.append([N, "ROOT", self._v(Nroots[I].E.real), self._v(Nroots[I].E.imag)])
            
    def _getRootsForN(self, N):
        for roots in self.allRoots:
            if roots.N == N:
                return roots
    
    def _v(self, num):
        return self._f(num, DISPLAY_DIFFPRESISION)
    
    def _f(self, num, precision):
        if abs(num) < pow(10,-precision):
            return "&lt;E-"+str(precision)
        else:
            if type(num) is str or type(num) is unicode:
                return num
            elif type(num) is QSType.mpmath.mpf:
                return QSType.mpmath.nstr(num, n=precision+QSType.mpIntDigits(num))
            else:
                return ('{:.'+str(precision)+'f}').format(num)
            
r = ResultsAnalyser(sys.argv[1],sys.argv[2],sys.argv[3])
r.createPoleTable()