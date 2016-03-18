from runBase import *
from Analytical.runArgsRange import *

seArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - All Signs", parents=[tcp_range])
args = seArgs.parse_args()

try:
  for i in [1.0,-1.0]:
    for j in [1.0,-1.0]:
      kCal = sm.kCalculator([args.t1_,args.t2_], 2.0, sm.K_SIGN, [i,j])
      sMats, ratSmat = getSmats(args, kCal, kCal)
except (DCException, sm.MatException) as inst:
  print str(inst)
  sys.exit()
