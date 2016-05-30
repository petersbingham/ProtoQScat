from runbase import *
from analytical.runargsrange import *

rtArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Range Test", parents=[tcp_range])
args = rtArgs.parse_args()

if QSMODE == MODE_MPMATH:
    print "Not supported for mpmath. Change mode in QSType."
else:
    try:
        kCal = sm.kCalculator([args.t1_,args.t2_], eneFactor=ENEFACTOR)
        anaSmat, ratSmat = getSmats(args, kCal, kCal)
        
        dEne = (args.eneEnd_-args.eneStart_) / float(args.eneSteps_)
        ene = args.eneStart_ + args.eneComplex_*1j 
        for i in range(0,args.eneSteps_+1,1):
            anaSmat.setEnergy(ene)
            ratSmat.setEnergy(ene)
            
            if np.allclose(ratSmat.getMatrix(), anaSmat.getMatrix(), num.MIN, num.MIN):
                print "\n" + str(ene) + " GOOD"
            else:
                print "\n" + str(ene) + " BAD"
            print str(ratSmat.getMatrix())
            print str(anaSmat.getMatrix())
            
            ene += dEne
        
    except (DCException, sm.MatException) as inst:
        print str(inst)
        sys.exit()