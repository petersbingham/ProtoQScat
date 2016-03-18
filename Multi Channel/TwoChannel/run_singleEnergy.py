from runBase import *
from Analytical.runArgsRange import *

seArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Single Energy", parents=[tcp_range])
seArgs.add_argument("ene_", help="Energy", type=complex)
args = seArgs.parse_args()

try:
  kCal = sm.kCalculator([args.t1_,args.t2_], 2.0)
  sMats, ratSmat = getSmats(args, kCal, kCal)
  ratSmat.setEnergy(args.ene_)  
  print str(ratSmat.getMatrix()) 
except (DCException, sm.MatException) as inst:
  print str(inst)
  sys.exit()
