from RatSMatWrap import *
import matplotlib.pyplot as plt
import argparse

parentArgs = argparse.ArgumentParser(description="Pyrazine Fit - Root find")
parentArgs.add_argument("subN_", help="Sub Set Number", type=int, nargs='?', default=None)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=None)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=None)
args = parentArgs.parse_args()

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  roots = kmat.findPolyRoots()
  for root in roots:
      print root
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()