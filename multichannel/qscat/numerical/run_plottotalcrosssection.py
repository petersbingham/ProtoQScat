from ratsmatwrap import *
from runargsplot import *

try:
    kmat = RatSMatWrap(FILENAME)
    anaxs = kmat.getTotalDiscreteXS()
    
    kmat = RatSMatWrap(FILENAME, args.subN_, args.subStart_, args.subEnd_)
    fitxs = kmat.getTotalDiscreteXS()
    fitxs.useMarker()
    ratTotxs = kmat.getTotalRatXS(eneStart=args.plotStart_, eneEnd=args.plotEnd_, eneComplexOffset=args.plotComplex_, eneSteps=args.steps_)
    
    anaxs.setDetails("R matrix", ['green'])
    fitxs.setDetails("Fit Points", ['red'])
    ratTotxs.setDetails("Fit Curve", ['blue'])
    
    sm.plotAll("Log", [anaxs, fitxs, ratTotxs], title="Fitted Numerical Data Cross Sections", xlabel="Electron Energy (rydbergs)", ylabel="Cross Section (bohr^2)") 
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()