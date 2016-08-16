import os
import sys
import datetime
import time
from sys import platform as _platform

BASEDIR = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,BASEDIR+'/../../../Utilities')
from general.file import *

RESULTSDIR = BASEDIR + sep() + "Results" + sep()
COEFFDIR = sep() + "CoefficientFiles" + sep()
ROOTSDIR = sep() + "Roots" + sep()
POLESDIR = sep() + "Poles"

LOGSDIR = BASEDIR + sep() + "Logs" + sep()

#Note to exceed the old max file path in win we must use \\?\ prefix, unicode and double slashes.
class ResultFileHandler:
    def __init__(self, sysName):
        self.sysPath = RESULTSDIR + sysName
        if not os.path.isdir(LOGSDIR):
            os.makedirs(LOGSDIR)
        self.logFile = open(LOGSDIR+str(datetime.datetime.now()).split('.')[0].replace(':','-')+".log", 'w+')
        self.startIndex = None
        self.endIndex = None
        self.numFits = None
        self.logTimes = {}
        
        #Following for the pole table. To keep track of max and min values
        self.numCmpStepsStart = None
        self.numCmpStepsEnd = None
        self.distFactorStart = None 
        self.distFactorEnd = None
    
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
        self.rootDir = base + ROOTSDIR
        self.rootPath = self.rootDir
        self._makeDir(self.rootPath)
        self.rootPath += self.getPostStr()
    
    def setPoleFindParameters(self, mode, numCmpSteps, distFactor, zeroVal):
        base = self._getRootsBaseString()
        self.polesDirName = POLESDIR + "_" + str(mode) + "_cfStep" + str(numCmpSteps) + "_dk" + str(distFactor) + "_zk" + str(zeroVal)
        self.polesDir = base + self.polesDirName + sep() 
        self.polesPath = self.polesDir
        self._makeDir(self.polesPath)
        self.polesPath += self.getPostStr()
        
        if self.numCmpStepsStart is None or numCmpSteps<self.numCmpStepsStart:
            self.numCmpStepsStart = numCmpSteps
        if self.numCmpStepsEnd is None or numCmpSteps>self.numCmpStepsEnd:
            self.numCmpStepsEnd = numCmpSteps
        if self.distFactorStart is None or distFactor<self.distFactorStart:
            self.distFactorStart = distFactor
        if self.distFactorEnd is None or distFactor>self.distFactorEnd:
            self.distFactorEnd = distFactor
        self.polesTableFile = base + sep() + str(mode) + "_cfStep" + str(self.numCmpStepsStart) + "-" + str(self.numCmpStepsEnd) + "_dk" + str(self.distFactorStart) + "-" + str(self.distFactorEnd) + "_zk" + str(zeroVal)

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
        self.writeLogStr(s)
        return fixPath(s)
    
    def getCoeffFileName(self, fit, ci, typeString, ext=".dat"):
        s = self.coeffPath + typeString + "_" + str(fit) + "_" + str(ci) + ext
        self.writeLogStr(s)
        return fixPath(s)

    def doesRootFileExist(self, ext=".dat"):
        return os.path.isfile(self.getRootFilePath(ext))

    def getRootDir(self):
        s = self.rootDir
        self.writeLogStr(s)
        return fixPath(s) 

    def getRootFilePath(self, ext=".dat"):
        s = self.rootPath + ext
        self.writeLogStr(s)
        return fixPath(s) 
    
    def getPoleDirName(self):
        s = self.polesDirName
        self.writeLogStr(s)
        return fixPath(s)

    def getPoleDir(self):
        s = self.polesDir
        self.writeLogStr(s)
        return fixPath(s)

    def getPoleFilePath(self, ext=".dat"):
        s = self.polesPath + ext
        self.writeLogStr(s)
        return fixPath(s)

    def getPoleTableParameters(self, ext=".tab"):
        s = self.polesTableFile + ext
        self.writeLogStr(s)
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
    
    def startLogAction(self, actionStr):
        self.writeLogStr("STARTING " + actionStr)
        self.logTimes[actionStr] = time.clock()
    
    def endLogAction(self, actionStr):
        self.writeLogStr("FINISHED " + actionStr + " Time Taken: " + str(time.clock()-self.logTimes[actionStr]))
    
    def writeLogStr(self, string):
        s = str(datetime.datetime.now().time()) + "\t" + string + "\n"
        self.logFile.write(s)
        self.logFile.flush()
    
    def _makeDir(self, path):
        if not os.path.isdir(fixPath(path)):
            os.makedirs(fixPath(path))