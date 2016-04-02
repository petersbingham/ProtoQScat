from RatSMatWrap import *
import argparse

parentArgs = argparse.ArgumentParser(description="Pyrazine Fit - Cross section")
parentArgs.add_argument("subN_", help="Sub Set Number", type=int, nargs='?', default=None)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=None)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=None)
args = parentArgs.parse_args()

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  totxs = kmat.getTotalDiscreteXS()
  ratTotxs = kmat.getTotalRatXS()
  
  totxs.setDetails("", ['red'])
  ratTotxs.setDetails("", ['blue'])
  
  sm.plotAll("Lin", [totxs, ratTotxs]) 
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()