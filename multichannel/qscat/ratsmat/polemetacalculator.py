import copy

from polefinder import *
from poleconverger import *
from general import *

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
            self._writePolePrevalenceTable(poleSetsDict, resultFileHandler)
            
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

    def _writePolePrevalenceTable(self, poleSetsDict, resultFileHandler):
        for cfstep in poleSetsDict:
            tabHeader = ["pole.E.real", "pole.E.imag", "Prevalence"]
            tabHeader.append(str(cfstep))
            poleSetsList = poleSetsDict[cfstep]
            uniquePoleSets = []
            totalContributingPoleTables = 0
            for poleSets in poleSetsList:
                if len(poleSets) > 0:
                    totalContributingPoleTables += 1
                    totalPoleCnt = reduce(lambda x,y: x+y, map(lambda poleSet: self._getNumPolesInPoleSet(poleSet), poleSets))
                    totalFactor = 0.0
                    for poleSet in poleSets:
                        i = self._getUniquePoleSetIndex(uniquePoleSets, poleSet)
                        newFactor = float(self._getNumPolesInPoleSet(poleSet))/totalPoleCnt
                        totalFactor += newFactor
                        if i == -1:
                            uniquePoleSets.append( [poleSet, newFactor] )
                        else:
                            oldFactor = uniquePoleSets[i][1]
                            uniquePoleSets[i] = [poleSet, oldFactor + newFactor] #Update pole set
                            #uniquePoleSets[i][1] = oldFactor + newFactor #Use origional pole set
            
            tabValues = []
            uniquePoleSets.sort(key=lambda x: x[1], reverse=True)       
            for uniquePoleSet in uniquePoleSets:
                Nmax = self._getMaxNInPoleSet(uniquePoleSet[0])
                prevalence = str(uniquePoleSet[1]/totalContributingPoleTables) + NOTABULATEFORMAT
                tabValues.append([formatRoot(uniquePoleSet[0][Nmax].E.real), formatRoot(uniquePoleSet[0][Nmax].E.imag), prevalence])
                
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
    
    def _getNumPolesInPoleSet(self, poleSet):
        return sum(not pole.isLost for pole in poleSet.values())
    
    def _getMaxNInPoleSet(self, poleSet):
        Nmax = 0
        for N in poleSet:
            pole = poleSet[N]
            if not pole.isLost and N>Nmax:
                Nmax = N
        return Nmax
                       