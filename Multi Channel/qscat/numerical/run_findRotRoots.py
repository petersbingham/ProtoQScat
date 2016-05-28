from RatSMatWrap import *
import matplotlib.pyplot as plt
from runBase import *

childArgs = argparse.ArgumentParser(description="Pyrazine Fit - All Sign Root find", parents=[parentArgs])
childArgs.add_argument("sene_", help="start energy", type=complex)
args = childArgs.parse_args()

FORWIKI = True

def getRootStr(root):
    if root is not None:
        return "{:.8f}".format(root)
    else:
        return "none\t"

print "\n\n"
first = True
try:
    kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
    roots = kmat.findConjRoots(args.sene_)
    print getRootStr(roots[0])+" || "+getRootStr(roots[1])+" ||"
                                
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()