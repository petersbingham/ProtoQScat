from runbase import *
from analytical.runargsrange import *

spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Pole find", parents=[tcp_range])
spArgs.add_argument("startE_", help="Start Energy", type=complex)
args = spArgs.parse_args()


try:
    kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, massMult=MASSMULT)
    anaSmat, ratSmat = getSmats(args, kCal, kCal)
  
    print "   " + str(ratSmat.findERoot(args.startE_)) + "   at +ve"
    print "   " + str(ratSmat.findERoot(args.startE_.real-1.0j*args.startE_.imag)) + "   at -ve"

except (DCException, sm.MatException) as inst:
    print str(inst)
    sys.exit()