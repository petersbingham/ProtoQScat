from ratsmatwrap import *
import matplotlib.pyplot as plt
import argparse
parentArgs = argparse.ArgumentParser(description="Numercal Data Fit Routine", add_help=False)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=-1)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=-1)
parentArgs.add_argument("sene_", help="start energy", type=complex)
args = parentArgs.parse_args()

NUM_FND = 0
          
def _findKnownRoot(N):
    kmat = RatSMatWrap(FILENAME,N,args.subStart_,args.subEnd_,kfitSigns=[1.0]*len(THRESHOLDS))
    starting = args.sene_
    print str(starting) + "   " + str(kmat.findERoot(starting)) + "\n\n"
        
try:
    for N in range (4,32,2):
        _findKnownRoot(N)
            
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()