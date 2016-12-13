from runargscommon import *

rootArgs = argparse.ArgumentParser(description="Two Channel Radial Well - Find Bound Roots", parents=[parentArgs])
args = rootArgs.parse_args()

kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0,1.0], massMult=MASSMULT)  
mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
smat = Smat(args.r0_, mats)
    
NUM_STEPS = 100

  
for i in range(1,NUM_STEPS):
    findStart = -1.0*args.v1_ + (args.v1_*float(i))/float(NUM_STEPS)
    try:
        root = str(smat.findRoot(findStart))
    except ValueError:
        root = "none" 
    print str(findStart) + "   " + str(root)
