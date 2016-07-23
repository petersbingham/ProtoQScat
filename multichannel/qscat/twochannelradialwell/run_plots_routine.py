from runargsplot import *
from runbase import *

IMAG = False

sm.setSubPlots(1,1,xlabel="Total Energy (hartrees)", ylabel="Imag" if IMAG else "Real")

kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, eneFactor=ENEFACTOR)
anaSmat = getAnaSmat(args, kCal)
ratSmat = getRatSmat(args, anaSmat, kCal, kCal)
doMatPlot(args, ratSmat, IMAG, "Approximated,  Imag Offset:" + str(args.plotComplex_), str(kCal))
