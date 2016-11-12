import copy

from polefinder import *
from poleconverger import *
from general import *
import general.numerical as num

class PoleMetaCalculator:
    def __init__(self, startIndex, endIndex, offset, mode, cfsteps, distFactors, zeroValExp, Nmin, Nmax, clusterSize, clusterLimit):
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
        self.clusterLimit = clusterLimit

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
            self._writePolePrevalenceTable(poleSetsDict, resultFileHandler)
            
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

    def _writePolePrevalenceTable(self, poleSetsDict, resultFileHandler):
        for cfstep in poleSetsDict:
            tabHeader = ["pole.E.real", "pole.E.imag", "Prevalence"]
            tabHeader.append(str(cfstep))
            poleSetsList = poleSetsDict[cfstep]
            uniquePoleLists = []
            for poleSets in poleSetsList:
                if len(poleSets) > 0:
                    totalPoleCnt = reduce(lambda x,y: x+y, map(lambda poleSet: self._getNumPolesInPoleSet(poleSet), poleSets))
                    for poleSet in poleSets:
                        i = self._getuniquePoleListIndex(uniquePoleLists, poleSet)
                        newFactor = float(self._getNumPolesInPoleSet(poleSet))/totalPoleCnt
                        if i == -1:
                            uniquePoleLists.append( (poleSet, newFactor) )
                        else:
                            oldFactor = uniquePoleLists[i][1]
                            uniquePoleLists[i] = (poleSet, oldFactor + newFactor)
            
            
            '''
            tabValues = []
            uniquePoleLists.sort(key=lambda x: x[1], reverse=True)       
            for uniquePoleList in uniquePoleLists:
                Nmax = self._getMaxNInPoleSet(uniquePoleList[0])
                prevalence = str(uniquePoleList[1]/len(poleSetsList)) + NOTABULATEFORMAT
                tabValues.append([formatRoot(uniquePoleList[0][Nmax].E.real), formatRoot(uniquePoleList[0][Nmax].E.imag), prevalence])
            outStr = getFormattedHTMLTable(tabValues, tabHeader, floatFmtFigs=DISPLAY_DIFFPRECISION, stralign="center", numalign="center", border=True)
            with open(resultFileHandler.getPolePrevalenceTablePath(cfstep, self.clusterSize, self.clusterLimit), 'w+') as f:
                f.write(outStr)
            continue
            '''
            
            
            poleValues = []     
            for uniquePoleList in uniquePoleLists:
                Nmax = self._getMaxNInPoleSet(uniquePoleList[0])
                prevalence = uniquePoleList[1]/len(poleSetsList)
                poleValues.append([uniquePoleList[0][Nmax], prevalence])

            ratCmp = num.RationalCompare(zeroValue=10**-self.zeroValExp, distFactor=self.clusterSize)                
            groupedPoleValuesList = []
            for poleValue in poleValues:
                found = False
                for groupedPoleValues in groupedPoleValuesList:
                    for groupedPoleValue in groupedPoleValues:
                        if ratCmp.isClose(poleValue[0].k, groupedPoleValue[0].k):
                            groupedPoleValues.append(poleValue)
                            found = True
                            break
                    if found:
                        break
                if not found:
                    groupedPoleValuesList.append([poleValue])
            
            totalPrevalence = 0.0
            poleValues = []    
            rejectedPoleValues = []    
            for groupedPoleValues in groupedPoleValuesList:
                for groupedPoleValue in groupedPoleValues:
                    if len(groupedPoleValues) <= self.clusterLimit:
                        poleValues.append(groupedPoleValue)
                        totalPrevalence += groupedPoleValue[1]
                    else:
                        rejectedPoleValues.append(groupedPoleValue)
       
            tabValues = []
            poleValues.sort(key=lambda x: x[1], reverse=True)
            for poleValues in poleValues:
                prevalence = str(poleValues[1]/totalPrevalence) + NOTABULATEFORMAT
                tabValues.append([formatRoot(poleValues[0].E.real), formatRoot(poleValues[0].E.imag), prevalence])
            outStr = getFormattedHTMLTable(tabValues, tabHeader, floatFmtFigs=DISPLAY_DIFFPRECISION, stralign="center", numalign="center", border=True)
            with open(resultFileHandler.getPolePrevalenceTablePath(cfstep, self.clusterSize, self.clusterLimit), 'w+') as f:
                f.write(outStr)

            rejectedTabValues = []
            for poleValues in rejectedPoleValues:
                rejectedTabValues.append([formatRoot(poleValues[0].E.real), formatRoot(poleValues[0].E.imag)])
            outStr = getFormattedHTMLTable(rejectedTabValues, tabHeader, floatFmtFigs=DISPLAY_DIFFPRECISION, stralign="center", numalign="center", border=True)
            with open(resultFileHandler.getRejectedPolePrevalenceTablePath(cfstep, self.clusterSize, self.clusterLimit), 'w+') as f:
                f.write(outStr)
    
    def _getuniquePoleListIndex(self, uniquePoleLists, poleSet):
        for i in range(len(uniquePoleLists)):
            uniquePoleList = self._getPolesInPoleSet(uniquePoleLists[i][0])
            found = True
            for N in sorted(poleSet.keys()):
                pole = poleSet[N]
                if pole.isLost:
                    break
                elif N not in uniquePoleList.keys() or pole not in uniquePoleList.values():
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
                       