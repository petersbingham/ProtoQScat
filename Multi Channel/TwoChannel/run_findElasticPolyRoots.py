from runBase import *
from Analytical.runArgsRange import *

spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Poly Root find", parents=[tcp_range])
args = spArgs.parse_args()

try:
    kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0,1.0], eneFactor=ENEFACTOR)
    getPolyRoots(args, kCal, kCal)
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()