import os
import sys
import datetime
import time
from sys import platform as _platform

basedir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,basedir+'/../../../Utilities')
from general.file import *

REPLACESTR = "REPLACETHIS"

BASEDIR = basedir

RESULTSDIR = BASEDIR + sep() + "Results" + sep()
COEFFDIR = sep() + "CoefficientFiles" + sep()
ROOTSDIR = sep() + "Roots" + sep()
ROOTSREJECTEDDIR = sep() + "Roots-Rejected" + sep()
POLESDIR = sep() + "Poles" + sep()
METACALSDIR = sep() + "MetaCals" + sep()

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
        self.distThresholdStart = None 
        self.distThresholdEnd = None
    
        self.cleanRootRoutine = None

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
        base = self._getBaseRootsPath()
        self.rootDir = base + ROOTSDIR
        self.rootPath = self.rootDir
        self._makeDir(self.rootPath)
        self.rootPath += self.getPostStr()
    
    def setCleanRootParameter(self, cleanRootRoutine):
        self.cleanRootRoutine = cleanRootRoutine
        if self.cleanRootRoutine is not None:
            self.cleanRootRoutine = str(cleanRootRoutine)
            base = self._getCleanRootsPath()
            
            self.cleanRootDir = base + ROOTSDIR
            self.cleanRootPath = self.cleanRootDir
            self._makeDir(self.cleanRootPath)
            self.cleanRootPath += self.getPostStr()
            
            self.rejectRootDir = base + ROOTSREJECTEDDIR
            self.rejectRootPath = self.rejectRootDir
            self._makeDir(self.rejectRootPath)
            self.rejectRootPath += self.getPostStr()
        
    def setPoleFindParameters(self, mode, numCmpSteps, distThreshold, zeroVal):
        self.mode = mode
        self.zeroVal = zeroVal
        
        base = self._getRootsPath()
        self.polesDirName = str(self.mode) + "_cfStep" + str(numCmpSteps) + "_dk" + str(distThreshold) + "_zk" + str(self.zeroVal)
        self.polesDir = base + POLESDIR + self.polesDirName + sep() 
        self.polesPath = self.polesDir
        self._makeDir(self.polesPath)
        self.polesPath += self.getPostStr()
        
        if self.numCmpStepsStart is None or numCmpSteps<self.numCmpStepsStart:
            self.numCmpStepsStart = numCmpSteps
        if self.numCmpStepsEnd is None or numCmpSteps>self.numCmpStepsEnd:
            self.numCmpStepsEnd = numCmpSteps
        if self.distThresholdStart is None or distThreshold<self.distThresholdStart:
            self.distThresholdStart = distThreshold
        if self.distThresholdEnd is None or distThreshold>self.distThresholdEnd:
            self.distThresholdEnd = distThreshold
    
    def setPoleMetaCalcParameters(self, amalgThreshold, Nmin, Nmax):        
        auxFilesRange = "_Nmin="+str(Nmin)+"_Nmax="+str(Nmax)
        auxFilesStrStart = str(self.mode) + auxFilesRange + "_cfStep"
        
        if self.distThresholdStart == self.distThresholdEnd:
            dkStr = "_dk" + str(self.distThresholdStart)
        else:
            dkStr = "_dk" + str(self.distThresholdStart) + "-" + str(self.distThresholdEnd)
        
        if self.numCmpStepsStart == self.numCmpStepsEnd:
            cfStepStr = auxFilesStrStart + str(self.numCmpStepsStart)
        else:
            cfStepStr = auxFilesStrStart + str(self.numCmpStepsStart) + "-" + str(self.numCmpStepsEnd)
        auxFilesStr1 = cfStepStr + dkStr + "_zk" + str(self.zeroVal)
        base = self._getRootsPath()
        self._makeDir(base + METACALSDIR)
        auxFilesStr2 = auxFilesStrStart + REPLACESTR + dkStr + "_ak" + str(amalgThreshold) + "_zk" + str(self.zeroVal)
        
        self.polesCountFile =  base + METACALSDIR + "COUNTS-" + auxFilesStr1
        self.polesPrevalenceFile =  base + METACALSDIR + auxFilesStr2 + "-PREVALENCE"
        self.combinedPolesPrevalenceFile =  base + METACALSDIR + auxFilesStr2 + "-PREVALENCECOMB"

    def getPostStr(self):
        if self.startIndex == None:
            self.startIndex = 0
        if self.endIndex == None:
            self.endIndex = self.fitSize-1
        return "N=" + str(self.fitSize) + "_" + "S=" + str(self.startIndex) + "_" + "E=" + str((self.endIndex+1*self.numFits)-1)
    
    def _getRootsInitPath(self):
        return self._getBaseString() + sep() + "ROOTS-" + self.rootFindRoutine + sep()
    
    def _getRootsPath(self):
        if self.cleanRootRoutine is None:
            return self._getBaseRootsPath()
        else:
            return self._getCleanRootsPath()
    
    def _getBaseRootsPath(self):
        return self._getRootsInitPath() + "COMPLETE"
    
    def _getCleanRootsPath(self):
        return self._getRootsInitPath() + "CLEAN-" + self.cleanRootRoutine
        
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

    def useCleanRoots(self):
        return self.cleanRootRoutine is not None

    def doesRootFileExist(self, ext=".dat"):
        return os.path.isfile(self.getRootFilePath(ext))

    def doesCleanRootFileExist(self, ext=".dat"):
        return self.useCleanRoots() and os.path.isfile(self.getCleanRootFilePath(ext))

    def getRootDir(self):
        s = self.rootDir
        self.writeLogStr(s)
        return fixPath(s) 

    def getRootFilePath(self, ext=".dat"):
        s = self.rootPath + ext
        self.writeLogStr(s)
        return fixPath(s) 

    def getCleanRootDir(self):
        s = self.cleanRootDir
        self.writeLogStr(s)
        return fixPath(s) 

    def getCleanRootFilePath(self, ext=".dat"):
        s = self.cleanRootPath + ext
        self.writeLogStr(s)
        return fixPath(s) 

    def getRejectRootFilePath(self, ext=".dat"):
        s = self.rejectRootPath + ext
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

    def getPoleCountTablePath(self, ext=".tab"):
        s = self.polesCountFile + ext
        self.writeLogStr(s)
        return fixPath(s)

    def getPolePrevalenceTablePath(self, distThreshold, ext=".tab"):
        s = self.polesPrevalenceFile.replace(REPLACESTR, str(distThreshold)) + ext
        self.writeLogStr(s)
        return fixPath(s)

    def getCombinedPolePrevalenceTablePath(self, distThreshold, ext=".tab"):
        s = self.combinedPolesPrevalenceFile.replace(REPLACESTR, str(distThreshold)) + ext
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