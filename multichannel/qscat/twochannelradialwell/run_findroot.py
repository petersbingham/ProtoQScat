from runbase import *
from analytical.runargsrange import *

spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Pole find", parents=[tcp_range])
spArgs.add_argument("startE_", help="Start Energy", type=complex)
args = spArgs.parse_args()

ANA_ROOT = 1.8315168862-0.0290733625j

#anakCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_ROT, massMult=MASSMULT)
#fitkCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, massMult=MASSMULT)

anakCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0,1.0], massMult=MASSMULT)
fitkCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0,1.0], massMult=MASSMULT)
ratkCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0,-1.0], massMult=MASSMULT)

anaSmat, ratSmat = getSmats(args, anakCal, fitkCal)
ratSmat.kCal = ratkCal
rootVal = ratSmat.findERoot(args.startE_)
print str(rootVal) + "   at +ve"

rootVal = ratSmat.findERoot(args.startE_.real - args.startE_.imag*1.0j)
print str(rootVal) + "   at -ve"

