import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')

from runbase import *

from analytical.runargsrange import *
spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Poly Root find", parents=[tcp_range])
spArgs.add_argument("mode_", help="Mode", type=int)
spArgs.add_argument("cfSteps_", help="Compare Steps", type=int)
spArgs.add_argument("distThreshold_", help="Distinguish Threshold", type=float)
spArgs.add_argument("zeroValExp_", help="Zero Value Precision", type=int)
spArgs.add_argument("cmpValue_", help="Compare Value", type=complex, nargs='?', default=None)

args = spArgs.parse_args()

MODE_DOUBLE = 0
MODE_INC = 1

def _doPoleFind(mode):
    try:
        kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0,1.0], massMult=MASSMULT)
        getPolyRoots(args, kCal, kCal, mode=mode, Nmax=30, cmpValue=args.cmpValue_)
    except (sm.MatException) as inst:
        raise inst
 
if args.mode_ == MODE_DOUBLE:
    _doPoleFind(DOUBLE_N)
elif args.mode_ == MODE_INC:
    _doPoleFind(INC_N)