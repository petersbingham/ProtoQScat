from ratsmatwrap import *
from runargsplot import *

try:
    kmat = RatSMatWrap(FILENAME)
    anaes = kmat.getDiscreteEigenSum()
    
    kmat = RatSMatWrap(FILENAME, args.subN_, args.subStart_, args.subEnd_, kfitSigns=[1.0]*kmat.numChannels)
    fites = kmat.getDiscreteEigenSum()
    fites.useMarker()
    ratTotes = kmat.getRatEigenSum(eneStart=args.plotStart_, eneEnd=args.plotEnd_, eneComplexOffset=args.offset_, eneSteps=args.steps_)
    
    anaes.setDetails("R-matrix", ['green'])
    fites.setDetails("Fit Points", ['red'])
    ratTotes.setDetails("Fit Curve", ['blue'])
    #, fites, ratTotes
    sm.plotSum("Lin", [anaes, fites, ratTotes], title="Fitted Numerical Data Eigenphase sums for N="+str(args.subN_), xlabel="Electron Energy (rydbergs)", ylabel="Total Phase") 
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()