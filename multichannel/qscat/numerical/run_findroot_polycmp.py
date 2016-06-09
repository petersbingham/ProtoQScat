from ratsmatwrap import *
import matplotlib.pyplot as plt
import argparse
parentArgs = argparse.ArgumentParser(description="Numercal Data Fit Routine", add_help=False)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=-1)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=-1)
args = parentArgs.parse_args()

KNOWN = True

NUM_FND = 0

def _findRoot(N, starting):
    kmat = RatSMatWrap(FILENAME,N,args.subStart_,args.subEnd_, kfitSigns=[1.0]*len(THRESHOLDS))
    new = round(starting.real,3) + round(starting.imag,4)*1.0j
    root = kmat.findRoot_(new)
    print str(new) + "   " + str(root) + "\n\n"
          
def _findKnownRoot(N):
    kmat = RatSMatWrap(FILENAME,N,args.subStart_,args.subEnd_, kfitSigns=[1.0]*len(THRESHOLDS))
    starting = KNOWN_POLE
    print str(starting) + "   " + str(kmat.findRoot(starting)) + "\n\n"
        
try:
    if KNOWN:
        for N in range (4,32,2):
            _findKnownRoot(N)
    else:
        i = 0
        for N in range (4,32,2):
            _findKnownRoot(N, STARTING_POSS[i])
            i += 1
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()