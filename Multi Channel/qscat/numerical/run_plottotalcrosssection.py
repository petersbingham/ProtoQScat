from RatSMatWrap import *
from runArgsPlot import *

try:
  kmat = RatSMatWrap("fort.19")
  anaxs = kmat.getTotalDiscreteXS()

  kmat = RatSMatWrap("fort.19", args.subN_, args.subStart_, args.subEnd_)
  fitxs = kmat.getTotalDiscreteXS()
  fitxs.useMarker()
  ratTotxs = kmat.getTotalRatXS(eneStart=args.plotStart_, eneEnd=args.plotEnd_, eneComplexOffset=args.plotComplex_, eneSteps=args.steps_)

  anaxs.setDetails("Analytical", ['green'])
  fitxs.setDetails("Fit Points", ['red'])
  ratTotxs.setDetails("Fit Curve", ['blue'])

  sm.plotAll("Lin", [anaxs, fitxs, ratTotxs], title="Fitted Pyrazine Cross Sections", xlabel="Electron Energy (rydbergs)", ylabel="Cross Section (bohr^2)") 
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()