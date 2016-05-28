from runArgsCommon import *

rootArgs = argparse.ArgumentParser(description="Two Channel Radial Well - Root Find", parents=[parentArgs])
rootArgs.add_argument("findStart_", help="Start Energy", type=complex)
args = rootArgs.parse_args()

for i in [1.0,-1.0]:
  for j in [1.0,-1.0]:
    kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=[i,j], eneFactor=ENEFACTOR)  
    mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
    smat = Smat(args.r0_, mats)
    try:
        root1 = str(smat.findRoot(args.findStart_))
    except ValueError:
        root1 = "none" 
    try:
        root2 = str(smat.findRoot(args.findStart_.real-1.0j*args.findStart_.imag))
    except ValueError:
        root2 = "none"
        
    print str([i,j]) + "\t" + root1 + "\t" + root2