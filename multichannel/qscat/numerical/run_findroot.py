from ratsmatwrap import *
import matplotlib.pyplot as plt
from runbase import *

childArgs = argparse.ArgumentParser(description="Numercal Data Fit - Root find", parents=[parentArgs])
childArgs.add_argument("sene_", help="start energy", type=complex)
args = childArgs.parse_args()

try:
    kmat = RatSMatWrap(FILENAME,args.subN_,args.subStart_,args.subEnd_)
    print kmat.findConjRoots(args.sene_, 1.0)
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()