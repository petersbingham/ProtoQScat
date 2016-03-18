import sys
sys.path.append("../Utilities")
import Scattering.Matrices as sm
import General.Numerical as num
from RatSMat import *

ZEROVALUE = 1E-7

class PoleFinder:
    def __init__(self, smats, kCal, resultsFolder, fitName, kConversionFactor, startIndex, endIndex, offset, distFactor, cmpValue=None):
        self.smats = smats
        self.kCal = kCal
        self.fitName = fitName
        self.kConversionFactor = kConversionFactor
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.offset = offset
        idStr = str(startIndex)+"_"+str(endIndex)+"_"+str(offset)+"_"+str(distFactor)+"_"+str(self.kCal)+".dat"
        self.file_coeff = open(resultsFolder+"/Output_Coeff_"+idStr, 'w')
        self.file_poles = open(resultsFolder+"/Output_"+idStr, 'w')
        self.ratCmp = num.RationalCompare(ZEROVALUE, distFactor)
        self.allPoles = []
        self.lastRoots = []
        self.cmpValue = cmpValue

        rootss = self._doubleM()
        self.file_coeff.close()
        self.file_poles.close()

    def _doubleM(self):
        Nmax = 64
        N = 4
        while N <= Nmax:
            roots = self._getMroots(N)
            self._locatePoles(roots)
            print str(len(roots)) + " roots."
            N = 2*N
        return roots

    def _getMroots(self, N):
        smats, step, endIndex = sm.decimate(self.smats, self.startIndex+self.offset, self.endIndex+self.offset, N)
        self.file_poles.write("\n")
        self._printSep2(self.file_poles)
        writeStr = "N=%d, Emin=%d, Emax=%d, step=%d, stepOff=%d\n" % (N,self.startIndex,endIndex,step,self.offset)
        self.file_poles.write(writeStr)
        ratSmat = RatSMat(smats, self.kCal.k, fitName=self.fitName)
        return ratSmat.findPolyRoots(self.kConversionFactor, False)
    
    def _calEnergy(self, k):
        return complex((1.0/self.kConversionFactor)*k**2)

    def _locatePoles(self, roots):
        newPoles = []
        numNewPoles = 0
        
        if self.cmpValue is not None:
            for j in range(len(roots)):
                eneRoot = self._calEnergy(roots[j])
                diff = abs(self.cmpValue-eneRoot)
                if j==0 or diff<closestDiff:
                    closestDiff = diff
                    closestIndex = j
            
        for i in range(len(roots)):
            root = roots[i]
            endStr1 = ""
            endStr2 = ""
            for j in range(len(self.lastRoots)):
              lastRoot = self.lastRoots[j]
              cdiff = self.ratCmp.getComplexDiff(root, lastRoot)
              if self.ratCmp.checkComplexDiff(cdiff):
                newPoles.append(root)
                numNewPoles += 1
                endStr1 = " diff[%d %d] = %.14f%+.14fi" % (i, j, cdiff.real, cdiff.imag)
            if self.cmpValue is not None and closestIndex==i:
                endStr2 = " @<****>@"
            eneRoot = self._calEnergy(root)
            writeStr = "Root_k[%d]=%.14f%+.14fi\tRoot_E[%d]=%.14f%+.14fi\t%s%s\n" % (i,root.real,root.imag,i,eneRoot.real,eneRoot.imag,endStr1,endStr2)
            self.file_poles.write(writeStr)
        self.lastRoots = roots

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

    def _printSep1(self, file):
        file.write("@<---------------------------------------------------------------------------------------------->@\n")

    def _printSep2(self, file):
        file.write("@<**********************************************************************************************>@\n")