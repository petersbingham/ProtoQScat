import os
import sys
import datetime
from sys import platform as _platform

BASEDIR = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,BASEDIR+'/../../../Utilities')
from general.file import *

RESULTSDIR = BASEDIR + sep() + "Results" + sep()
COEFFDIR = sep() + "CoefficientFiles" + sep()
ROOTSDIR = sep() + "Roots" + sep()
POLESDIR = sep() + "Poles"

#Note to exceed the old max file path in win we must use \\?\ prefix, unicode and double slashes.
class ResultFileHandler:
    def __init__(self, sysName):
        self.sysPath = RESULTSDIR + sysName
        self.logFile = open(BASEDIR+"/path.log", 'w')
        self.startIndex = None
        self.endIndex = None
        self.numFits = None
    
    def setFitInfo(self, numFits, fitSize):
        self.numFits = numFits
        self.fitSize = fitSize
        self.numCoeffs = fitSize/2 + 1
    
    def setDeciminationInfo(self, startIndex, endIndex):
        self.startIndex = startIndex
        self.endIndex = endIndex
    
    def setCoeffRoutine(self, coeffRoutine):
        self.coeffRoutinePath = self.sysPath + sep() + "COEFFS-" + coeffRoutine
        
        base = self._getBaseString()
        self.coeffPath = base+COEFFDIR+self.getPostStr()+sep()
        self._makeDir(self.coeffPath)
    
    def setRootFindRoutine(self, rootFindRoutine):
        self.rootFindRoutine = rootFindRoutine
        base = self._getRootsBaseString()
        self.rootPath = base + ROOTSDIR
        self._makeDir(self.rootPath)
        self.rootPath += self.getPostStr()
    
    def setPoleFindParameters(self, mode, numCmpSteps, distFactor, zeroVal):
        base = self._getRootsBaseString()
        self.polesPath = base + POLESDIR + "_" + str(mode) + "_cfStep" + str(numCmpSteps) + "_dk" + str(distFactor) + "_zk" + str(zeroVal) + sep() 
        self._makeDir(self.polesPath)
        self.polesPath += self.getPostStr()
    
    def getPostStr(self):
        if self.startIndex == None:
            self.startIndex = 0
        if self.endIndex == None:
            self.endIndex = self.fitSize-1
        return "N=" + str(self.fitSize) + "_" + "S=" + str(self.startIndex) + "_" + "E=" + str((self.endIndex+1*self.numFits)-1)
    
    def _getRootsBaseString(self):
        return self._getBaseString() + sep() + "ROOTS-" + self.rootFindRoutine
    
    def _getBaseString(self):
        if self.numFits == 1:
            return self.coeffRoutinePath + sep() + "SingleFit"
        else:
            return self.coeffRoutinePath + sep() + str(self.numFits) + "Fits"
                     
    def doCoeffFilesExist(self):
        if os.path.isdir(self.coeffPath):
            for file in self._getAfitNames() + self._getBfitNames():
                if not os.path.isfile(file):
                    return False
            return True
        return False  
    
    def getCoeffFilePath(self):
        s = self.coeffPath
        self._printLogStr(s)
        return fixPath(s)
    
    def getCoeffFileName(self, fit, ci, typeString, ext=".dat"):
        s = self.coeffPath + typeString + "_" + str(fit) + "_" + str(ci) + ext
        self._printLogStr(s)
        return fixPath(s)

    def doesRootFileExist(self, ext=".dat"):
        return os.path.isfile(self.getRootFilePath(ext))

    def getRootFilePath(self, ext=".dat"):
        s = self.rootPath + ext
        self._printLogStr(s)
        return fixPath(s) 

    def getPoleFilePath(self, ext=".dat"):
        s = self.polesPath + ext
        self._printLogStr(s)
        return fixPath(s)

    def _getAfitNames(self): 
        return self._getFitNames("A")
      
    def _getBfitNames(self): 
        return self._getFitNames("B")
    
    def _getFitNames(self, typeString): 
        fitNames = []
        for fit in range(self.numFits):
            for ci in range(0, self.numCoeffs):
                fitNames.append(self.getCoeffFileName(fit, ci, typeString))
        return fitNames  
    
    def _printLogStr(self, string):
        s = str(datetime.datetime.now().time()) + "\t" + string + "\n"
        self.logFile.write(s)
    
    def _makeDir(self, path):
        if not os.path.isdir(fixPath(path)):
            os.makedirs(fixPath(path))