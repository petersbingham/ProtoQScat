from ratsmatwrap import *
import matplotlib.pyplot as plt
from runbase import *

childArgs = argparse.ArgumentParser(description="Numercal Data Fit - All Sign Root find", parents=[parentArgs])
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
    kmat = RatSMatWrap(FILENAME,args.subN_,args.subStart_,args.subEnd_)
    roots = kmat.findERoot_Conj(args.sene_)
    print getRootStr(roots[0])+" || "+getRootStr(roots[1])+" ||"
                                
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()