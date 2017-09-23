from ratsmatwrap import *
import matplotlib.pyplot as plt
from runbase import *
import general.numerical as num

childArgs = argparse.ArgumentParser(description="Numercal Data Fit - All Sign Root find", parents=[parentArgs])
childArgs.add_argument("sene_", help="start energy", type=complex)
args = childArgs.parse_args()

FORWIKI = True

def getRootStr(root):
    if root is not None:
        return "{:.8f}".format(root)
    else:
        return "none\t"

print "\n\n"
first = True
try:
    kperms = num.getPermutations([1.0,-1.0], len(THRESHOLDS))
    for kfitSigns in kperms:
        if FORWIKI:
            print "\nFit Signs:" + str(kfitSigns)
        else:
            print "\n************* Fit Signs:" + str(kfitSigns) + "*************"
        if first:
            first = False
            if FORWIKI:
                print "||~ ||~ 1st Starting Point ||~ 2nd Starting Point ||"
                print "||~ Jost Signs ||~ " + str(args.sene_) + " ||~ " + str(args.sene_.real-1.0j*args.sene_.imag)+" ||"
            else:
                print "\t\t\t1st Starting Point\t\t2nd Starting Point"
                print "Jost Signs\t\t" + str(args.sene_) + "\t\t\t" + str(args.sene_.real-1.0j*args.sene_.imag)
        for ksigns in kperms:
            kmat = RatSMatWrap(FILENAME,args.subN_,args.subStart_,args.subEnd_,kfitSigns,ksigns,True)
            roots = kmat.findERoot_Conj(args.sene_)
            if FORWIKI:
                print "|| "+str(ksigns)+" || "+getRootStr(roots[0])+" || "+getRootStr(roots[1])+" ||"
            else:
                print str(ksigns)+"   \t"+getRootStr(roots[0])+"   \t"+getRootStr(roots[1])
                                
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()