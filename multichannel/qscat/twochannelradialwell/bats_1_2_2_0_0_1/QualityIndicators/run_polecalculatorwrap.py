import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../..')
sys.path.insert(0,base+'/../../..')
sys.path.insert(0,base+'/../../../ratsmat')
sys.path.insert(0,os.getcwd()) #We assume that the specific kreader and description (below) will be here.

from resultfilehandler import *
from polemetacalculator import *
from runbase import *
from analytical.runargsrange import *

MODE_ALLSIGNS_DOUBLE = 0
MODE_ALLSIGNS_INC = 1
MODE_ROT_DOUBLE = 2
MODE_ROT_INC = 3
MODE_POS_DOUBLE = 4
MODE_POS_INC = 5
MODE_COMP_DOUBLE = 6
MODE_COMP_INC = 7

import argparse
parentArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - QI calculator", parents=[tcp_range])
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

def _doPoleFind(kCal, mode):
    calculateQIs(args, kCal, kCal, mode)

def _polesForAllSigns(mode):
    kperms = num.getPermutations([1.0,-1.0], 2)
    for kperm in kperms:
        kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=kperm, massMult=MASSMULT)
        _doPoleFind(kCal, mode)

def _polesForRot(mode):
    kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_ROT, massMult=MASSMULT)
    _doPoleFind(kCal, mode)

def _polesForPos(mode):
    kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0]*2, massMult=MASSMULT)
    _doPoleFind(kCal, mode)

def _polesForComp(mode):
    kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, massMult=MASSMULT)
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
elif args.mode_ == MODE_COMP_DOUBLE:
    _polesForComp(DOUBLE_N)
elif args.mode_ == MODE_COMP_INC:
    _polesForComp(INC_N)