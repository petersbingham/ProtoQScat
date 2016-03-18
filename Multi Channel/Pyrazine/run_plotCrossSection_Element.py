from RatSMatWrap import *
import matplotlib.pyplot as plt
import argparse

parentArgs = argparse.ArgumentParser(description="Pyrazine Fit - Cross Section Plot per element")
parentArgs.add_argument("m_", help="m element", type=int)
parentArgs.add_argument("n_", help="n element", type=int)
parentArgs.add_argument("subN_", help="Sub Set Number", type=int, nargs='?', default=None)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=None)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=None)
args = parentArgs.parse_args()

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  xs = kmat.getDiscreteXS()
  ratxs = kmat.getRatXS()
  
  xs.setDetails("Raw", ['red'])
  ratxs.setDetails("Ana", ['blue'])
  
  sm.plotSingle("Lin", [xs, ratxs], args.m_, args.n_, False, "Pyrazine Cross Sections", 'Electron Energy (rydbergs)', 'Cross Section (bohr^2)')
  
  plt.show()
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()