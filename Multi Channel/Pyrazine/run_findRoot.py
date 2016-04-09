from RatSMatWrap import *
import matplotlib.pyplot as plt
from runBase import *

childArgs = argparse.ArgumentParser(description="Pyrazine Fit - Root find", parents=[parentArgs])
childArgs.add_argument("sene_", help="start energy", type=complex)
args = childArgs.parse_args()

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  print kmat.findConjRoots(args.sene_, 1.0)
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()