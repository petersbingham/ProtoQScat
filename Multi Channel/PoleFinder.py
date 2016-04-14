import sys
import os
sys.path.append("../Utilities")
import Scattering.Matrices as sm
import General.Numerical as num
from RatSMat import *

ZEROVALUE = 1E-7
DOUBLE_N = 0
INC_N = 1

class PoleFinder:
    def __init__(self, smats, kCal, resultsFolder, fitName, kConversionFactor, startIndex, endIndex, offset, distFactor, numCmpSteps=1, cmpValue=None, mode=DOUBLE_N):
        self.smats = smats
        self.kCal = kCal
        self.fitName = fitName
        self.kConversionFactor = kConversionFactor
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.offset = offset
        self.numCmpSteps = numCmpSteps
        idStr = str(startIndex)+"_"+str(endIndex)+"_"+str(offset)+"_"+str(distFactor)+"_"+str(self.kCal)+".dat"
        self.file_coeff = open(resultsFolder+"/Output_Coeff_"+idStr, 'w')
        self.file_poles = open(resultsFolder+"/Output_"+idStr, 'w')
        self.ratCmp = num.RationalCompare(ZEROVALUE, distFactor)
        self.allPoles = []
        self.allRoots = []
        self.allNs = []
        self.cmpValue = cmpValue

        if mode == DOUBLE_N:
            self._doubleN()
        else:
            self._incN()
        
        self.file_coeff.close()
        self.file_poles.close()

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
        roots = self._getNroots(N)
        self._locatePoles(roots)
        self.allNs.append(N)
        print str(len(roots)) + " roots.\n\n"

    def _getNroots(self, N):
        actualStartIndex = self.startIndex+self.offset
        smats, step, actualEndIndex, startEne, endEne = sm.decimate(self.smats, actualStartIndex, self.endIndex+self.offset, N)
        self.file_poles.write("\n")
        self._printSep2(self.file_poles)
        decStr = "N=%d, Emin=%d(%f), Emax=%d(%f), step=%d" % (N,actualStartIndex,startEne,actualEndIndex,endEne,step)
        print "Decimation:"
        print "  "+decStr
        self.file_poles.write(decStr+"\n")
        ratSmat = RatSMat(smats, self.kCal, fitName=self.fitName)
        return ratSmat.findPolyRoots(self.kConversionFactor, False)

    def _locatePoles(self, roots):
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
        numNewPoles = 0
        if len(self.allRoots) >= self.numCmpSteps:  
            for i in range(len(roots)):
                root = roots[i]
                isPole = True
                endStr = ""
                for k in range(len(self.allRoots)-self.numCmpSteps, len(self.allRoots)): #Look at the last sets
                    cmpRootSet = self.allRoots[k]
                    for j in range(len(cmpRootSet)):
                        cmpRoot = cmpRootSet[j]
                        cdiff = self.ratCmp.getComplexDiff(root, cmpRoot)
                        if self.ratCmp.checkComplexDiff(cdiff):
                            endStr += " with N=%d[%d]: diff = %.14f%+.14fi" % (self.allNs[k], j, cdiff.real, cdiff.imag)
                            break
                        if j==len(cmpRootSet)-1:
                            endStr = ""
                            isPole = False
                    if not isPole:
                        break
                if isPole:        
                    newPoles.append(root)
                    numNewPoles += 1
                self._printRoot(i, root, endStr, closestIndex)
        else:
            for i in range(len(roots)):
                self._printRoot(i, roots[i], "", closestIndex)

        self.allRoots.append(roots)
        
        self._printSep1(self.file_poles)

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
        for newPole in newPoles:
          isNew = True
          for i in range(len(self.allPoles)):
            pole = self.allPoles[i]
            if self.ratCmp.isClose(newPole, pole):
              isNew = False
              break

          if isNew:
            self.allPoles.append(newPole)
            #Record when we start adding new poles
            if newIndex == -1:
              newIndex = len(self.allPoles)-1
          else:
            #If it's not new then just update the value
            self.allPoles[i] = newPole

        for i in range(len(self.allPoles)):
          endStr = ""
          enePole = self._calEnergy(self.allPoles[i])
          if newIndex!=-1 and i>=newIndex:
            endStr = "NEW"
            
          for j in range(len(lostIndices)):
            if i == lostIndices[j]:
              endStr = "LOST"
              break

          writeStr = "Pole_k[%d]=%.14f%+.14fi\tPole_E[%d]=%.14f%+.14fi\t%s\n" % (i,self.allPoles[i].real,self.allPoles[i].imag,i,enePole.real,enePole.imag,endStr)
          self.file_poles.write(writeStr)
    
    def _calEnergy(self, k):
        return complex((1.0/self.kConversionFactor)*k**2)

    def _printRoot(self, i, root, endStr, closestIndex):
        endStr2 = ""
        if self.cmpValue is not None and closestIndex==i:
            endStr2 = " @<****>@"
        eneRoot = self._calEnergy(root)
        writeStr = "Root_k[%d]=%.14f%+.14fi\tRoot_E[%d]=%.14f%+.14fi\t%s%s\n" % (i,root.real,root.imag,i,eneRoot.real,eneRoot.imag,endStr,endStr2)
        self.file_poles.write(writeStr)
                
    def _printSep1(self, file):
        file.write("@<---------------------------------------------------------------------------------------------->@\n")

    def _printSep2(self, file):
        file.write("@<**********************************************************************************************>@\n")