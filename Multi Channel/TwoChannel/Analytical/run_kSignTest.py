from runArgsCommon import *

seArgs = argparse.ArgumentParser(description="Two Channel Radial Well - Single Energy", parents=[parentArgs])
args = seArgs.parse_args()

ene = 1.9
print "Energy = " + str(ene)
kCal = sm.kCalculator([args.t1_,args.t2_], 2.0, sm.K_SIGN, [1.0,1.0])
mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
smat = Smat(args.r0_, mats)  
smat.setEnergy(ene)
print "\nS:"
print str(smat)

kCal = sm.kCalculator([args.t1_,args.t2_], 2.0, sm.K_SIGN, [-1.0,-1.0])
mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
smat = Smat(args.r0_, mats)  
smat.setEnergy(ene)
print "\nS:"
print str(smat)
