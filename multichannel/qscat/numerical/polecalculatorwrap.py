import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')
sys.path.insert(0,os.getcwd()) #We assume that the specific kreader and description (below) will be here.

from matreader import *
from sysdesc import *
from polemetacalculator import *

CLUSTERSIZE = 10**-3
CLUSTERLMT = 2

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
parentArgs.add_argument("distFactors_", help="Distinguish Factor", type=str)
parentArgs.add_argument("zeroValExp_", help="Zero Value Precision", type=int)
parentArgs.add_argument("Nmin_", help="Starting N value", type=int)
parentArgs.add_argument("Nmax_", help="Ending N value", type=int)
parentArgs.add_argument("cmpIndex_", help="Compare Index", type=int, nargs='?', default=0)
args = parentArgs.parse_args()

from scriptparameters import *

def _doPoleFind(kCal, mode):
    kmats = readkMats(FILENAME)
    smats = sm.getSfromKmatrices(kmats,NUMCHANNELS)
    if args.endIndex_ == -1:
        args.endIndex_ = len(smats)-1
    resultFileHandler = getFileHandler(kCal, args.startIndex_, args.endIndex_)
    
    cfsteps = map(int, args.cfSteps_.split(','))
    distFactors = map(float, args.distFactors_.split(','))
    p = PoleMetaCalculator(args.startIndex_, args.endIndex_, args.offset_, mode, cfsteps, distFactors, args.zeroValExp_, args.Nmin_, args.Nmax_, CLUSTERSIZE, CLUSTERLMT)
    cmpPole = RMATRIX_POLES[args.cmpIndex_] if args.cmpIndex_<len(RMATRIX_POLES) else None
    p.doPoleCalculations(smats, resultFileHandler, kCal, mode, cmpPole)

def _polesForAllSigns(mode):
    kperms = num.getPermutations([1.0,-1.0], len(THRESHOLDS))
    for kperm in kperms:
        kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=kperm, eneFactor=ENEFACTOR)
        _doPoleFind(kCal, mode)

def _polesForRot(mode):
    kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_ROT, eneFactor=ENEFACTOR)
    _doPoleFind(kCal, mode)

def _polesForPos(mode):
    kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=[1.0]*len(THRESHOLDS), eneFactor=ENEFACTOR)
    _doPoleFind(kCal, mode)

if args.mode_ == MODE_ALLSIGNS_DOUBLE:
    _polesForAllSigns(DOUBLE_N)
elif args.mode_ == MODE_ALLSIGNS_INC:
    _polesForAllSigns(INC_N)
elif args.mode_ == MODE_ROT_DOUBLE:
    _polesForRot(DOUBLE_N)
elif args.mode_ == MODE_ROT_INC:
    _polesForRot(INC_N)
elif args.mode_ == MODE_POS_DOUBLE:
    _polesForPos(DOUBLE_N)
elif args.mode_ == MODE_POS_INC:
    _polesForPos(INC_N)