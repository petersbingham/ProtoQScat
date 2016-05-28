from runBase import *
from Analytical.runArgsRange import *

spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Pole find", parents=[tcp_range])
spArgs.add_argument("startE_", help="Start Energy", type=complex)
args = spArgs.parse_args()


for i in [1.0,-1.0]:
  for j in [1.0,-1.0]:
    try:
        kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[i,j], eneFactor=ENEFACTOR)
        anaSmat, ratSmat = getSmats(args, kCal, kCal)
      
        try:
            print str([i,j]) + "   " + str(ratSmat.findRoot(args.startE_)) + "   at +ve"
        except ValueError:
            try:
                print str([i,j]) + "   " + str(ratSmat.findRoot(-1.0*args.startE_)) + "   at -ve"
            except ValueError:
                print str([i,j]) + "   Val Err"
      
    except (DCException, sm.MatException) as inst:
      print str(inst)
      sys.exit()