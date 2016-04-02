from RatSMatWrap import *
import matplotlib.pyplot as plt
from runBase import *

childArgs = argparse.ArgumentParser(description="Pyrazine Fit - All Sign Root find", parents=[parentArgs])
childArgs.add_argument("sene_", help="start energy", type=complex)
args = childArgs.parse_args()

def getRootStr(root):
    if root is not None:
        return "{:.8f}".format(root)
    else:
        return "none\t"

print "\n\n"
print "\t\t\t1st Starting Point\t\t2nd Starting Point"
print "Jost Signs\t\t" + str(args.sene_) + "\t\t\t" + str(args.sene_.real-1.0j*args.sene_.imag)
try:
    for i in [1.0,-1.0]:
        for j in [1.0,-1.0]:
            for k in [1.0,-1.0]:
                kfitSigns = [i,j,k]
                print "\n************* Fit Signs:" + str(kfitSigns) + "*************"
                for ii in [1.0,-1.0]:
                    for jj in [1.0,-1.0]:
                        for kk in [1.0,-1.0]:
                            ksigns = [ii,jj,kk]
                            kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_,kfitSigns,ksigns,True)
                            roots = kmat.findConjRoots(args.sene_)
                            print str(ksigns)+"   \t"+getRootStr(roots[0])+"   \t"+getRootStr(roots[1])
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()