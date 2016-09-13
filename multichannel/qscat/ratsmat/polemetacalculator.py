import copy

from polefinder import *
from poleconverger import *
from general import *

class PoleMetaCalculator:
    def __init__(self, startIndex, endIndex, offset, mode, cfsteps, distFactors, zeroValExp, Nmin, Nmax):
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.offset = offset
        self.mode = mode
        self.cfsteps = cfsteps
        self.distFactors = distFactors
        self.zeroValExp = zeroValExp
        self.Nmin = Nmin
        self.Nmax = Nmax

    def doPoleCalculations(self, smats, resultFileHandler, kCal, mode, cmpPole=None):
        tabList = []
        for distFactor in self.distFactors:
            for cfStep in self.cfsteps:
                p = PoleFinder(copy.deepcopy(smats), kCal, resultFileHandler, self.startIndex, self.endIndex, self.offset, distFactor, cfStep, cmpPole, mode, zeroValExp=self.zeroValExp, Nmin=self.Nmin, Nmax=self.Nmax)
                tabList.append((p.NmaxTotPoles, p.NmaxLostPoles))
                r = PoleConverger(resultFileHandler)
                r.createPoleTable()
        self._writePoleCountTables(tabList, resultFileHandler)
            
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