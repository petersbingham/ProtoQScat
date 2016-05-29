from RatSMatWrap import *
from runbase import *
args = parentArgs.parse_args()

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  xs = kmat.getDiscreteXS()
  ratxs = kmat.getRatXS()
  
  xs.setDetails("", ['red','green','blue','purple','orange','brown','cyan','pink','grey'])
  ratxs.setDetails("", ['red','green','blue','purple','orange','brown','cyan','pink','grey'])
  
  sm.plotAll("Lin", [xs], False, "Analytical")
  sm.plotAll("Lin", [ratxs], False, "Approximated")
  
  plt.show()
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()