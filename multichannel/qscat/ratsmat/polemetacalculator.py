import copy
import math

from polefinder import *
from poleconverger import *
from general import *

CALCULATIONS = ["Q2", "Q3", "Q4", "Q5", "Q6", "Q7"]


class PoleMetaCalculator:
    def __init__(self, startIndex, endIndex, offset, mode, cfsteps, distFactors, relaxFactor, zeroValExp, Nmin, Nmax):
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.offset = offset
        self.mode = mode
        self.cfsteps = cfsteps
        self.distFactors = distFactors
        self.relaxFactor = relaxFactor
        self.zeroValExp = zeroValExp
        self.Nmin = Nmin
        self.Nmax = Nmax

    def doPoleCalculations(self, smats, resultFileHandler, kCal, mode, cmpPole=None):
        if len(smats) <= self.endIndex:
            raise Exception("Specified End Index outside range.")
        if self.startIndex < 0:
            raise Exception("Specified Start Index less than zero.")
        tabList = []
        poleSetsDict = {}
        self.errState = False
        for cfStep in self.cfsteps:
            poleSetsDict[cfStep] = []
            for distFactor in sorted(self.distFactors, reverse=True):  #Want sorted for the prevalence, since each pole across the N in each pole set should be a subset of the same pole for a higer distFactor 
                pf = PoleFinder(copy.deepcopy(smats), kCal, resultFileHandler, self.startIndex, self.endIndex, self.offset, distFactor, self.relaxFactor, cfStep, cmpPole, mode, zeroValExp=self.zeroValExp, Nmin=self.Nmin, Nmax=self.Nmax)
                self.errState = self.errState | pf.errState
                if not self.errState:
                    tabList.append((pf.NmaxTotPoles, pf.NmaxLostPoles))
                    pc = PoleConverger(resultFileHandler, self.Nmin, self.Nmax)
                    pc.createPoleTable()
                    poleSetsDict[cfStep].append(pc.poleSets)
        if not self.errState:
            self._writePoleCountTables(tabList, resultFileHandler)
            self._writePoleCalculationTable(poleSetsDict, resultFileHandler)
            
    def _writePoleCountTables(self, tabList, resultFileHandler):
        tabHeader = ["dk"]
        for cfstep in self.cfsteps:
            if cfstep == 1:
                tabHeader.append(str(cfstep)+" Step")
            else:
                tabHeader.append(str(cfstep)+" Steps")
        
        tabValues = []
        for id in range(len(self.distFactors)):
            tabRow = ["{:.2E}".format(self.distFactors[id])]
            for ic in range(len(self.cfsteps)):
                index = ic * len(self.distFactors) + id
                if tabList[index][1] > 0:
                    tabRow.append(str(tabList[index][0])+"("+str(tabList[index][1])+")")
                else:
                    tabRow.append(str(tabList[index][0]))
            tabValues.append(tabRow)
                
        outStr = getFormattedHTMLTable(tabValues, tabHeader, stralign="center", numalign="center", border=True)
        with open(resultFileHandler.getPoleCountTablePath(), 'w+') as f:
            f.write(outStr)
        print outStr

    def _writePoleCalculationTable(self, poleSetsDict, resultFileHandler):   
        for cfstep in poleSetsDict:
            poleSetsList = poleSetsDict[cfstep]
            uniquePoleSets = []
            kTOT = 0.0
            lenPiSumkSumi = 0.0
            totPoleCnts = 0.0
            for poleSets in poleSetsList:
                if len(poleSets) > 0:
                    kTOT += 1
                    lenPiSumk = reduce(lambda x,y: x+y, map(lambda poleSet: self._getLenpi(poleSet), poleSets))
                    lenPiSumkSumi += lenPiSumk
                    for poleSet in poleSets:
                        i = self._getUniquePoleSetIndex(uniquePoleSets, poleSet)
                        
                        lenPi = self._getLenpi(poleSet)
                        q1_inter = float(lenPi)/lenPiSumk
                        q2_inter = lenPi
                        q5_inter = 1.0
                        totPoleCnts += 1.0
                        
                        if i == -1:
                            uniquePoleSets.append( [poleSet, q1_inter, q2_inter, q5_inter] )
                        else:
                            q1_inter_old = uniquePoleSets[i][1]
                            q2_inter_old = uniquePoleSets[i][2]
                            q5_inter_old = uniquePoleSets[i][3]
                            uniquePoleSets[i] = [poleSet, q1_inter_old+q1_inter, q2_inter_old+q2_inter, q5_inter_old+q5_inter] #Update pole set
            self._writeTable(resultFileHandler, cfstep, uniquePoleSets, kTOT, lenPiSumkSumi, totPoleCnts)
            
    def _writeTable(self, resultFileHandler, cfstep, uniquePoleSets, kTOT, lenPiSumkSumi, totPoleCnts):
        if self.mode == DOUBLE_N:
            nTOT = math.log(self.Nmax,2.0) - math.log(self.Nmin,2.0) + 1.0
        else:
            nTOT = self.Nmax/2.0 - self.Nmin/2.0 + 1.0        
        tabValues = []
        uniquePoleSets.sort(key=lambda x: x[1], reverse=True)
        
        tabHeader = ["pole.E.real", "pole.E.imag", "Total n Range", "Number ks Located"]
        #tabHeader = ["pole.E.real", "pole.E.imag", "Q1", "Q2", "Q3", "Q4", "Q5", "Q6", "Q7"]
        for uniquePoleSet in uniquePoleSets:
            Nmax = self._getMaxNInPoleSet(uniquePoleSet[0])
            q1 = str(uniquePoleSet[1]/kTOT) + NOTABULATEFORMAT
            q2 = str(uniquePoleSet[2]) + NOTABULATEFORMAT
            q3 = str(uniquePoleSet[2]/nTOT) + NOTABULATEFORMAT
            q4 = str(uniquePoleSet[2]/lenPiSumkSumi) + NOTABULATEFORMAT
            q5 = str(uniquePoleSet[3]) + NOTABULATEFORMAT
            q6 = str(uniquePoleSet[3]/nTOT) + NOTABULATEFORMAT
            q7 = str(uniquePoleSet[3]/totPoleCnts) + NOTABULATEFORMAT
            tabValues.append([formatRoot(uniquePoleSet[0][Nmax].E.real), formatRoot(uniquePoleSet[0][Nmax].E.imag), q2, q5])
            #tabValues.append([formatRoot(uniquePoleSet[0][Nmax].E.real), formatRoot(uniquePoleSet[0][Nmax].E.imag), q1, q2, q3, q4, q5, q6, q7])
            
        outStr = getFormattedHTMLTable(tabValues, tabHeader, floatFmtFigs=DISPLAY_DIFFPRECISION, stralign="center", numalign="center", border=True)
        with open(resultFileHandler.getPolePrevalenceTablePath(cfstep), 'w+') as f:
            f.write(outStr)
    
    def _getUniquePoleSetIndex(self, uniquePoleSets, poleSet):
        for i in range(len(uniquePoleSets)):
            uniquePoleSet = self._getPolesInPoleSet(uniquePoleSets[i][0])
            found = True
            for N in sorted(poleSet.keys()):
                pole = poleSet[N]
                if pole.isLost:
                    break
                elif N not in uniquePoleSet.keys() or pole not in uniquePoleSet.values():
                    found = False
                    break
            if found:
                return i
        return -1
    
    def _getPolesInPoleSet(self, poleSet):
        nonLostPoles = {}
        for N in poleSet:
            pole = poleSet[N]
            if not pole.isLost:
                nonLostPoles[N] = pole
        return nonLostPoles
    
    def _getLenpi(self, poleSet):
        return sum(not pole.isLost for pole in poleSet.values())
    
    def _getMaxNInPoleSet(self, poleSet):
        Nmax = 0
        for N in poleSet:
            pole = poleSet[N]
            if not pole.isLost and N>Nmax:
                Nmax = N
        return Nmax
                       