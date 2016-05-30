from runargscommon import *

seArgs = argparse.ArgumentParser(description="Two Channel Radial Well - Single Energy", parents=[parentArgs])
seArgs.add_argument("ene_", help="Energy", type=complex)
args = seArgs.parse_args()

kCal = sm.kCalculator([args.t1_,args.t2_], eneFactor=ENEFACTOR)
mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
smat = Smat(args.r0_, mats)  
smat.setEnergy(args.ene_)
print "\nS:"
print str(smat)
print "\n|S|:"
print str(smat.abs())
print "\nS'S:"
print str(smat.uniOp())
if LIN_ALGEBRA:
    print "\na:"
    print str(mats.a)
    print "\na^2:"
    print str(mats.aSq)
print "\nA:"
print str(mats.A)
print "\nV:"
print str(mats.V)
print "\nK:"
print str(mats.K)
