import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')
sys.path.insert(0,os.getcwd()) #We assume that the specific kreader and description (below) will be here.

from matreader import *
from sysdesc import *
from polemetacalculator import *

MODE_ALLSIGNS_DOUBLE = 0
MODE_ALLSIGNS_INC = 1
MODE_ROT_DOUBLE = 2
MODE_ROT_INC = 3
MODE_POS_DOUBLE = 4
MODE_POS_INC = 5
MODE_COMP_DOUBLE = 6
MODE_COMP_INC = 7
MODE_NEG_DOUBLE = 8
MODE_NEG_INC = 9
MODE_CUST_DOUBLE = 10
MODE_CUST_INC = 11

import argparse
parentArgs = argparse.ArgumentParser(description="Numercal Data Fit - Pole find")
parentArgs.add_argument("startIndex_", help="Start Index", type=int)
parentArgs.add_argument("endIndex_", help="End Index", type=int)
parentArgs.add_argument("offset_", help="Offset", type=int)
parentArgs.add_argument("mode_", help="Mode", type=int)
parentArgs.add_argument("cfSteps_", help="Compare Steps", type=str)
parentArgs.add_argument("startingDistThreshold_", help="Starting Distinguish Threshold", type=float)
parentArgs.add_argument("amalgThreshold_", help="Amalgamation Threshold", type=float)
parentArgs.add_argument("zeroValExp_", help="Zero Value Precision", type=int)
parentArgs.add_argument("Nmin_", help="Starting N value", type=int)
parentArgs.add_argument("Nmax_", help="Ending N value", type=int)
parentArgs.add_argument("cmpIndex_", help="Compare Index", type=int, nargs='?', default=0)
args = parentArgs.parse_args()

from scriptparameters import *

kCAL=[1.0]*len(THRESHOLDS)

def _doPoleFind(ratkCal, mode):
    kmats = readkMats(FILENAME)
    smats = sm.getSfromKmatrices(kmats,NUMCHANNELS)
    if args.endIndex_ == -1:
        args.endIndex_ = len(smats)-1
    resultFileHandler = getFileHandler(kCAL, args.startIndex_, args.endIndex_)
    resultFileHandler.setRatsMatCalcStr(str(ratkCal))
    cfSteps = map(int, args.cfSteps_.split(','))
    p = PoleMetaCalculator(args.startIndex_, args.endIndex_, args.offset_, mode, cfSteps, args.startingDistThreshold_, args.amalgThreshold_, args.zeroValExp_, args.Nmin_, args.Nmax_, resultFileHandler)
    cmpPole = RMATRIX_POLES[args.cmpIndex_] if args.cmpIndex_<len(RMATRIX_POLES) else None
    
    p.doPoleCalculations(smats, kCAL, mode, ratkCal=ratkCal, cmpPole=cmpPole)

def _polesForAllSigns(mode):
    kperms = num.getPermutations([1.0,-1.0], len(THRESHOLDS))
    for kperm in kperms:
        ratkCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=kperm, massMult=MASSMULT)
        _doPoleFind(ratkCal, mode)

def _polesForRot(mode):
    ratkCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_ROT, massMult=MASSMULT)
    _doPoleFind(ratkCal, mode)

def _polesForPos(mode):
    ratkCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=[1.0]*len(THRESHOLDS), massMult=MASSMULT)
    _doPoleFind(ratkCal, mode)

def _polesForNeg(mode):
    ratkCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=[-1.0]*len(THRESHOLDS), massMult=MASSMULT)
    _doPoleFind(ratkCal, mode)

def _polesForCust(mode):
    ratkCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=[-1.0,1.0,-1.0], massMult=MASSMULT)
    _doPoleFind(ratkCal, mode)

def _polesForComp(mode):
    ratkCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_COMP, massMult=MASSMULT)
    _doPoleFind(ratkCal, mode)

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
elif args.mode_ == MODE_COMP_DOUBLE:
    _polesForComp(DOUBLE_N)
elif args.mode_ == MODE_COMP_INC:
    _polesForComp(INC_N)
elif args.mode_ == MODE_NEG_DOUBLE:
    _polesForNeg(DOUBLE_N)
elif args.mode_ == MODE_NEG_INC:
    _polesForNeg(INC_N)
elif args.mode_ == MODE_CUST_DOUBLE:
    _polesForCust(DOUBLE_N)
elif args.mode_ == MODE_CUST_INC:
    _polesForCust(INC_N)