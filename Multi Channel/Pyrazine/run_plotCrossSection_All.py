from RatSMatWrap import *
import argparse

parentArgs = argparse.ArgumentParser(description="Pyrazine Fit - Cross section")
parentArgs.add_argument("subN_", help="Sub Set Number", type=int, nargs='?', default=None)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=None)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=None)
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