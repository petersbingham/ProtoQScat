from runArgsCommon import *

rootArgs = argparse.ArgumentParser(description="Two Channel Radial Well - Root Find", parents=[parentArgs])
rootArgs.add_argument("findStart_", help="Start Energy", type=complex)
args = rootArgs.parse_args()

for i in [1.0,-1.0]:
  for j in [1.0,-1.0]:
    kCal = sm.kCalculator([args.t1_,args.t2_], 2.0, sm.K_SIGN, [i,j])  
    mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
    smat = Smat(args.r0_, mats)
    try:
        print str([i,j]) + "   " + str(smat.findRoot(args.findStart_))
    except ValueError:
        try:
            print str([i,j]) + "   " + str(smat.findRoot(-1.0*args.findStart_))
        except ValueError:
            print str([i,j]) + "   Val Err"