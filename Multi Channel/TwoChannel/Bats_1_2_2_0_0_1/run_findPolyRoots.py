import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')

from runBase import *

from Analytical.runArgsRange import *
spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Poly Root find", parents=[tcp_range])
spArgs.add_argument("mode_", help="Mode", type=int)
spArgs.add_argument("cmpValue_", help="Compare Value", type=complex, nargs='?', default=None)
args = spArgs.parse_args()

MODE_DOUBLE = 0
MODE_INC = 1

def _doPoleFind(mode, dirName):
    path = "./Results/"+dirName
    if not os.path.exists(path):
        os.makedirs(path)
    try:
        kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0,1.0], eneFactor=ENEFACTOR)
        getPolyRoots(args, kCal, kCal, path, cmpValue=args.cmpValue_, mode=mode)
    except (sm.MatException) as inst:
      print str(inst)
      sys.exit()
 
if args.mode_ == MODE_DOUBLE:
    _doPoleFind(DOUBLE_N, "Double")
elif args.mode_ == MODE_INC:
    _doPoleFind(INC_N, "Inc")