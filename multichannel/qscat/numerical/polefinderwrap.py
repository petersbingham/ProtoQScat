import sys
import os
import copy
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')
sys.path.insert(0,os.getcwd()) #We assume that the specific kreader and description (below) will be here.

from matreader import *
from sysdesc import *
from polefinder import *
from resultsanalyser import *
from general import *

MODE_ALLSIGNS_DOUBLE = 0
MODE_ALLSIGNS_INC = 1
MODE_ROT_DOUBLE = 2
MODE_ROT_INC = 3
MODE_POS_DOUBLE = 4
MODE_POS_INC = 5

import argparse
parentArgs = argparse.ArgumentParser(description="Numercal Data Fit - Pole find")
parentArgs.add_argument("startIndex_", help="Start Index", type=int)
parentArgs.add_argument("endIndex_", help="End Index", type=int)
parentArgs.add_argument("offset_", help="Offset", type=int)
parentArgs.add_argument("mode_", help="Mode", type=int)
parentArgs.add_argument("cfSteps_", help="Compare Steps", type=str)
parentArgs.add_argument("distFactor_", help="Distinguish Factor", type=str)
parentArgs.add_argument("zeroValExp_", help="Zero Value Precision", type=int)
parentArgs.add_argument("Nmin_", help="Starting N value", type=int)
parentArgs.add_argument("Nmax_", help="Ending N value", type=int)
parentArgs.add_argument("cmpIndex_", help="Compare Index", type=int, nargs='?', default=0)
args = parentArgs.parse_args()

from scriptparameters import *

def _doPoleFind(kCal, mode, dirName):
    kmats = readkMats(FILENAME)
    smats = sm.getSfromKmatrices(kmats,NUMCHANNELS)
    resultFileHandler = getFileHandler(kCal, args.startIndex_, args.endIndex_)
    tabList = []
    for distFactor in map(float, args.distFactor_.split(',')):
        for cfstep in map(int, args.cfSteps_.split(',')):    
            cmpPole = RMATRIX_POLES[args.cmpIndex_] if args.cmpIndex_<len(RMATRIX_POLES) else None
            p = PoleFinder(copy.deepcopy(smats), kCal, resultFileHandler, args.startIndex_, args.endIndex_, args.offset_, distFactor, cfstep, cmpPole, mode, Nmin=args.Nmin_, Nmax=args.Nmax_)
            tabList.append((p.NmaxTotPoles, p.NmaxLostPoles))
            r = ResultsAnalyser(resultFileHandler)
            r.createPoleTable()
    _writePoleTables(tabList, resultFileHandler)
        
def _writePoleTables(tabList, resultFileHandler):
    splitcfSteps = args.cfSteps_.split(',')
    splitDistFactor = args.distFactor_.split(',')
    
    tabHeader = ["dk"]
    for cfstep in map(int, splitcfSteps):
        if cfstep == 1:
            tabHeader.append(str(cfstep)+" Step")
        else:
            tabHeader.append(str(cfstep)+" Steps")
    
    tabValues = []
    i=0
    for distFactor in map(float, splitDistFactor):
        tabRow = ["{:.2E}".format(distFactor)]
        for cfstep in map(int, splitcfSteps):
            if tabList[i][1] > 0:
                tabRow.append(str(tabList[i][0])+"("+str(tabList[i][1])+")")
            else:
                tabRow.append(str(tabList[i][0]))
            i += 1
        tabValues.append(tabRow)
            
    outStr = getFormattedHTMLTable(tabValues, tabHeader, stralign="center", numalign="center", border=True)
    with open(resultFileHandler.getPoleTableParameters(), 'w+') as f:
        f.write(outStr)
    print outStr

def _polesForAllSigns(mode, dirName):
    kperms = num.getPermutations([1.0,-1.0], len(THRESHOLDS))
    for kperm in kperms:
        kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=kperm, eneFactor=ENEFACTOR)
        _doPoleFind(kCal, mode, dirName)

def _polesForRot(mode, dirName):
    kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_ROT, eneFactor=ENEFACTOR)
    _doPoleFind(kCal, mode, dirName)

def _polesForPos(mode, dirName):
    kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=[1.0]*len(THRESHOLDS), eneFactor=ENEFACTOR)
    _doPoleFind(kCal, mode, dirName)

if args.mode_ == MODE_ALLSIGNS_DOUBLE:
    _polesForAllSigns(DOUBLE_N, "AllSigns Double")
elif args.mode_ == MODE_ALLSIGNS_INC:
    _polesForAllSigns(INC_N, "AllSigns Inc")
elif args.mode_ == MODE_ROT_DOUBLE:
    _polesForRot(DOUBLE_N, "Rot Double")
elif args.mode_ == MODE_ROT_INC:
    _polesForRot(INC_N, "Rot Inc")
elif args.mode_ == MODE_POS_DOUBLE:
    _polesForPos(DOUBLE_N, "Pos Double")
elif args.mode_ == MODE_POS_INC:
    _polesForPos(INC_N, "Pos Inc")