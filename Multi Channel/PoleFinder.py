import sys
import os
import traceback
sys.path.append("../Utilities")
import Scattering.Matrices as sm
import General.Numerical as num
from RatSMat import *
from General.QSType import *


ZEROVALEXP = 7
DOUBLE_N = 'doubleN'
INC_N = 'incN'
COMPLETE_STR = "Complete"

class Decimator():
    def __init__(self, startIndex, endIndex, offset, resultFileHandler):
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.offset = offset
        self.resultFileHandler = resultFileHandler
                 
    def decimate(self, sMats, N):
        actualStartIndex = self.startIndex+self.offset
        sMats, step, actualEndIndex, startEne, endEne = sm.decimate(sMats, actualStartIndex, self.endIndex+self.offset, N)
        self.resultFileHandler.setDeciminationInfo(actualStartIndex, actualEndIndex)
        try:
            decStr = "N=%d, Emin=%d(%f), Emax=%d(%f), step=%d" % (N,actualStartIndex,startEne,actualEndIndex,endEne,step)
        except TypeError:
            decStr = "N=%d, Emin=%d(%f+%fi), Emax=%d(%f+%fi), step=%d" % (N,actualStartIndex,startEne.real,startEne.imag,actualEndIndex,endEne.real,endEne.imag,step)
        print "Decimation:"
        print "  "+decStr
        return sMats, decStr     

class PoleFinder:
    def __init__(self, sMats, kCal, resultFileHandler, kConversionFactor, startIndex, endIndex, offset, distFactor, numCmpSteps=1, cmpValue=None, mode=DOUBLE_N, populateSmatCB=None, zeroValExp=ZEROVALEXP):
        self.sMats = sMats
        self.kCal = kCal
        self.resultFileHandler = resultFileHandler
        self.kConversionFactor = kConversionFactor
        self.decimator = Decimator(startIndex, endIndex, offset, resultFileHandler)
        self.numCmpSteps = numCmpSteps
        self.zeroValExp = zeroValExp
        self.zeroValue = 10**(-zeroValExp)
        self.ratCmp = num.RationalCompare(self.zeroValue, distFactor)
        self.allPoles = []
        self.allPolesInfoStrs = []
        self.allRoots = []
        self.allNs = []
        self.cmpValue = cmpValue
        self.populateSmatCB = populateSmatCB
        self.file_roots = None
        self.file_poles = None
        self.distFactor = distFactor
        self.mode = mode
        self.first = True
        
        if mode == DOUBLE_N:
            self._doubleN()
        else:
            self._incN()

    def _doubleN(self):
        Nmax = 64
        N = 4
        while N <= Nmax:
            self._doForN(N)
            N = 2*N

    def _incN(self):
        Nmax = 64
        N = 4
        while N <= Nmax:
            self._doForN(N)
            N = N+2

    def _doForN(self, N):
        try:
            if not self.first:
                print "\n"
            self.first = False
            roots = self._getNroots(N)
            self._locatePoles(roots, N)
            self.allNs.append(N)
        except Exception as inst:
            self.allRoots.append([])
            string = "Unhandled Exception: " + str(inst) + "\n"
            traceback.print_exc(file=sys.stdout)
            if self.file_roots is not None:
                self.file_roots.write(string)
            if self.file_poles is not None:
                self.file_poles.write(string)
        if self.file_roots is not None:
            self.file_roots.close()
            self.file_roots = None
        if self.file_poles is not None:
            self.file_poles.close()
            self.file_poles = None

    def _getNroots(self, N):
        sMats, decStr = self.decimator.decimate(self.sMats, N)
        ratSmat = RatSMat(sMats, self.kCal, resultFileHandler=self.resultFileHandler, doCalc=False)
        self.resultFileHandler.setPoleFindParameters(self.mode, self.numCmpSteps, self.distFactor, self.zeroValExp)
        roots = None
        if self.resultFileHandler.doesRootFileExist():
            ratSmat.coeffSolve.printCalStr(True)
            try:
                roots = self._readNroots(N)
                ratSmat.polyRootSolve.printCalStr(True)
                print "Loaded Roots for N=" + str(N) + ":"
            except Exception as e:
                print "Error reading roots will attempt to recalculate"
        if roots is None:
            if self.populateSmatCB is not None:
                self.populateSmatCB(sMats, self.sMats)
            ratSmat.doCalc()
            self.file_roots = open(self.resultFileHandler.getRootFilePath(), 'w')
            self.file_roots.write(decStr+"\n")
            roots = ratSmat.findPolyRoots(self.kConversionFactor, False)
            print "Calculated Roots for N=" + str(N) + ":"
        print "  " + str(len(roots)) + " roots."
        self.file_poles = open(self.resultFileHandler.getPoleFilePath(), 'w')
        self.file_poles.write(decStr+"\n")
        return roots

    def _readNroots(self, N):
        roots = []
        fileComplete = False
        with open(self.resultFileHandler.getRootFilePath(), 'r') as f:
            firstLine = True
            for line in f:        
                if COMPLETE_STR in line:
                    return roots
                elif not firstLine:
                    str = line[line.find('=')+1:line.find('i')]+'j'
                    roots.append(complex(str))
                firstLine = False
        raise Exception("Incomplete Root File")

    def _locatePoles(self, roots, N):
        #This is when we know the position of a pole and want to mark the closest root to this value in the output file.
        closestIndex = None
        if self.cmpValue is not None:
            for j in range(len(roots)):
                eneRoot = self._calEnergy(roots[j])
                diff = abs(self.cmpValue-eneRoot)
                if j==0 or diff<closestDiff:
                    closestDiff = diff
                    closestIndex = j

        newPoles = []
        newPolesInfoStrs = []
        numNewPoles = 0
        if len(self.allRoots) >= self.numCmpSteps:  
            for i in range(len(roots)):
                root = roots[i]
                isPole = True
                infoStr = "with N=%d[%d]" % (N, i)
                for k in range(len(self.allRoots)-self.numCmpSteps, len(self.allRoots)): #Look at the last sets
                    cmpRootSet = self.allRoots[k]
                    for j in range(len(cmpRootSet)):
                        cmpRoot = cmpRootSet[j]
                        cdiff = self.ratCmp.getComplexDiff(root, cmpRoot)
                        if self.ratCmp.checkComplexDiff(cdiff):
                            infoStr += " N=%d[%d]" % (self.allNs[k], j)
                            break
                        if j==len(cmpRootSet)-1:
                            isPole = False
                    if not isPole:
                        break
                if isPole:        
                    newPoles.append(root)
                    newPolesInfoStrs.append(infoStr)
                    numNewPoles += 1
                self._printRoot(i, root, closestIndex)
        else:
            for i in range(len(roots)):
                self._printRoot(i, roots[i], closestIndex)
        if self.file_roots is not None:
            self.file_roots.write(COMPLETE_STR)
        

        self.allRoots.append(roots)

        #Determine if lost pole. If so then note.
        lostIndices = []
        for i in range(len(self.allPoles)):
          pole = self.allPoles[i]
          lostPole = True;
          for newPole in newPoles:
            if self.ratCmp.isClose(pole, newPole):
              lostPole = False
              break
          if lostPole:
            lostIndices.append(i)

        #Determine if new pole. If so then add to self.allPoles and note.
        newIndex = -1
        for i in range(len(newPoles)):
          newPole = newPoles[i]
          newPolesInfoStr = newPolesInfoStrs[i]
          isNew = True
          for i in range(len(self.allPoles)):
            pole = self.allPoles[i]
            if self.ratCmp.isClose(newPole, pole):
              isNew = False
              break

          if isNew:
            self.allPoles.append(newPole)
            self.allPolesInfoStrs.append(newPolesInfoStr)
            #Record when we start adding new poles
            if newIndex == -1:
              newIndex = len(self.allPoles)-1
          else:
            #If it's not new then just update the value
            self.allPoles[i] = newPole
            self.allPolesInfoStrs[i] = newPolesInfoStr

        poles = 0
        newPoles = 0
        lostPoles = 0
        for i in range(len(self.allPoles)):
          poles += 1
          endStr = ""
          enePole = self._calEnergy(self.allPoles[i])
          if newIndex!=-1 and i>=newIndex:
            endStr = "NEW "
            newPoles += 1
            
          for j in range(len(lostIndices)):
            if i == lostIndices[j]:
              endStr = "LOST "
              poles -= 1
              lostPoles += 1
              break

          writeStr = ("Pole_k[%d]="+self._getComplexFormat()+"\tPole_E[%d]="+self._getComplexFormat()+"    \t%s%s\n") % (i,self.allPoles[i].real,self.allPoles[i].imag,i,enePole.real,enePole.imag,endStr,self.allPolesInfoStrs[i])
          self.file_poles.write(writeStr)
          
        self.file_poles.write(COMPLETE_STR)
        print "Poles calculated in mode " + self.mode + ", using df=" + str(self.distFactor)
        print "Calculated Poles for N=" + str(N) + ":"
        print "  " + str(poles) + " poles, of which " + str(newPoles) + " are new. " + str(lostPoles) + (" has" if lostPoles==1 else " have") + " been lost."
    
    def _calEnergy(self, k):
        return complex((1.0/self.kConversionFactor)*k**2)

    def _printRoot(self, i, root, closestIndex):
        endStr = ""
        if self.cmpValue is not None and closestIndex==i:
            endStr = " @<****>@"
        eneRoot = self._calEnergy(root)
        writeStr = ("Root_k[%d]="+self._getComplexFormat()+"\tRoot_E[%d]="+self._getComplexFormat()+"\t%s\n") % (i,root.real,root.imag,i,eneRoot.real,eneRoot.imag,endStr)
        if self.file_roots is not None:
            self.file_roots.write(writeStr)
                
    def _printSep1(self, file):
        file.write("@<---------------------------------------------------------------------------------------------->@\n")

    def _printSep2(self, file):
        file.write("@<**********************************************************************************************>@\n")
        
    def _getComplexFormat(self):
        if QSMODE == MODE_NORM:
            return "%.14f%+.14fi"
        else:
            return "%."+str(DPS)+"f%+."+str(DPS)+"fi"   