from RatSMatWrap import *
from runBase import *
args = parentArgs.parse_args()

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  totxs = kmat.getTotalDiscreteXS()
  totxs.useMarker()
  ratTotxs = kmat.getTotalRatXS()
  
  totxs.setDetails("", ['red'])
  ratTotxs.setDetails("", ['blue'])
  
  sm.plotAll("Lin", [totxs, ratTotxs]) 
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()