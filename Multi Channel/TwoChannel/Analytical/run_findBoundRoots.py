from runArgsCommon import *

rootArgs = argparse.ArgumentParser(description="Two Channel Radial Well - Find Bound Roots", parents=[parentArgs])
args = rootArgs.parse_args()

kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[1.0,1.0], eneFactor=ENEFACTOR)  
mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
smat = Smat(args.r0_, mats)
    
for i in range(1,100):
    findStart = -args.v1_ + float(i)/(args.v1_*100.0)
    try:
        root = str(smat.findRoot(findStart))
    except ValueError:
        root = "none" 
    print root