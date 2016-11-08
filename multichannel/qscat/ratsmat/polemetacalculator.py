import copy

from polefinder import *
from poleconverger import *
from general import *
import general.numerical as num

class PoleMetaCalculator:
    def __init__(self, startIndex, endIndex, offset, mode, cfsteps, distFactors, zeroValExp, Nmin, Nmax, clusterSize):
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.offset = offset
        self.mode = mode
        self.cfsteps = cfsteps
        self.distFactors = distFactors
        self.zeroValExp = zeroValExp
        self.Nmin = Nmin
        self.Nmax = Nmax
        self.clusterSize = clusterSize

    def doPoleCalculations(self, smats, resultFileHandler, kCal, mode, cmpPole=None):
        if len(smats) <= self.endIndex:
            raise Exception("Specified End Index outside range.")
        if self.startIndex < 0:
            raise Exception("Specified Start Index less than zero.")
        tabList = []
        poleSetsDict = {}
        self.errState = False
        pc = None
        for cfStep in self.cfsteps:
            poleSetsDict[cfStep] = []
            for distFactor in sorted(self.distFactors, reverse=True):  #Want sorted for the prevalence, since each pole across the N in each pole set should be a subset of the same pole for a higer distFactor 
                pf = PoleFinder(copy.deepcopy(smats), kCal, resultFileHandler, self.startIndex, self.endIndex, self.offset, distFactor, cfStep, cmpPole, mode, zeroValExp=self.zeroValExp, Nmin=self.Nmin, Nmax=self.Nmax)
                self.errState = self.errState | pf.errState
                if not self.errState:
                    tabList.append((pf.NmaxTotPoles, pf.NmaxLostPoles))
                    pc = PoleConverger(resultFileHandler, self.Nmin, self.Nmax)
                    pc.createPoleTable()
                    poleSetsDict[cfStep].append(pc.poleSets)
        if not self.errState and pc is not None:
            self._writePoleCountTables(tabList, resultFileHandler)
            self._writePolePrevalenceTable(pc.allRoots, poleSetsDict, resultFileHandler)
            
    def _writePoleCountTables(self, tabList, resultFileHandler):
        tabHeader = ["dk"]
        for cfstep in self.cfsteps:
            if cfstep == 1:
                tabHeader.append(str(cfstep)+" Step")
            else:
                tabHeader.append(str(cfstep)+" Steps")
        
        tabValues = []
        i=0
        for distFactor in self.distFactors:
            tabRow = ["{:.2E}".format(distFactor)]
            for cfstep in self.cfsteps:
                if tabList[i][1] > 0:
                    tabRow.append(str(tabList[i][0])+"("+str(tabList[i][1])+")")
                else:
                    tabRow.append(str(tabList[i][0]))
                i += 1
            tabValues.append(tabRow)
                
        outStr = getFormattedHTMLTable(tabValues, tabHeader, stralign="center", numalign="center", border=True)
        with open(resultFileHandler.getPoleCountTablePath(), 'w+') as f:
            f.write(outStr)
        print outStr

    def _writePolePrevalenceTable(self, allRoots, poleSetsDict, resultFileHandler):
        for cfstep in poleSetsDict:
            tabHeader = ["pole.E.real", "pole.E.imag", "Prevalence", "Cluster Factor", "Goodness"]
            tabHeader.append(str(cfstep))
            poleSetsList = poleSetsDict[cfstep]
            uniquePoleSets = []
            for poleSets in poleSetsList:
                if len(poleSets) > 0:
                    totalPoleCnt = reduce(lambda x,y: x+y, map(lambda poleSet: self._getNumPolesInPoleSet(poleSet), poleSets))
                    poleSetFactors = []
                    totalModFactor = 0.0
                    for poleSet in poleSets:
                        factor = float(self._getNumPolesInPoleSet(poleSet))/totalPoleCnt
                        clusterFac = self._getCumulativeClusterFactor(allRoots, poleSet)
                        modFactor = factor * clusterFac
                        totalModFactor += modFactor
                        poleSetFactors.append((poleSet, factor, clusterFac, modFactor))
                    
                    normFactor = 1.0 / totalModFactor
                    for poleSetFactor in poleSetFactors:
                        poleSet = poleSetFactor[0] 
                        factor = poleSetFactor[1]
                        clusterFac = poleSetFactor[2]
                        finalModFactor = poleSetFactor[3] * normFactor
                        i = self._getUniquePoleSetIndex(uniquePoleSets, poleSet)
                        if i == -1:
                            uniquePoleSets.append( (poleSet, factor, clusterFac, finalModFactor, 1) )
                        else:
                            cumulatingFactor = uniquePoleSets[i][1]
                            cumulatingClusterFac = uniquePoleSets[i][2]
                            cumulatingModFactor = uniquePoleSets[i][3]
                            cnt = uniquePoleSets[i][4]
                            uniquePoleSets[i] = (poleSet, cumulatingFactor+factor, (cumulatingClusterFac*cnt+clusterFac)/(cnt+1), cumulatingModFactor+finalModFactor, cnt+1)
            
            tabValues = []
            uniquePoleSets.sort(key=lambda x: x[3], reverse=True)       
            for uniquePoleSet in uniquePoleSets:
                Nmax = self._getMaxNInPoleSet(uniquePoleSet[0])
                prevalence = str(uniquePoleSet[1]/len(poleSetsList)) + NOTABULATEFORMAT
                clusterFactor = str(uniquePoleSet[2]) + NOTABULATEFORMAT
                modPrevalence = str(uniquePoleSet[3]/len(poleSetsList)) + NOTABULATEFORMAT
                tabValues.append([formatRoot(uniquePoleSet[0][Nmax].E.real), formatRoot(uniquePoleSet[0][Nmax].E.imag), prevalence, clusterFactor, modPrevalence])
                
            outStr = getFormattedHTMLTable(tabValues, tabHeader, floatFmtFigs=DISPLAY_DIFFPRECISION, stralign="center", numalign="center", border=True)
            with open(resultFileHandler.getPolePrevalenceTablePath(cfstep, self.clusterSize), 'w+') as f:
                f.write(outStr)
                
    def _getCumulativeClusterFactor(self, allRoots, poleSet):
        zeroValue = 10**-self.zeroValExp
        ratCmp = num.RationalCompare(zeroValue, self.clusterSize)
        clusterFact = 1.0
        for N in poleSet:
            clusterCnt = 0.0
            if not poleSet[N].isLost:
                roots = self._getRoots(allRoots, N)
                for root in roots:
                    if ratCmp.isClose(root.k, poleSet[N].k):
                        clusterCnt += 1.0
                if clusterCnt == 0.0:
                    raise Exception("Zero cluster count!") #Should never be here since we know the pole will always be here at the least.
                clusterFact = clusterFact / clusterCnt
        return clusterFact
    
    def _getHighestNInPolesSet(self, poleSet):
        highN = 0
        for N in poleSet:
            if not poleSet[N].isLost:
                if N > highN:
                    highN = N
        return highN
    
    def _getRoots(self, allRoots, N):
        for roots in allRoots:
            if roots.N == N:
                return roots
        raise Exception("Could not find roots for pole set!") #Should never be here
    
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
                       