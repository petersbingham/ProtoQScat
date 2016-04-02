from RatSMatWrap import *
import matplotlib.pyplot as plt
import argparse

parentArgs = argparse.ArgumentParser(description="Pyrazine Fit - Cross section")
parentArgs.add_argument("subN_", help="Sub Set Number", type=int, nargs='?', default=None)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=None)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=None)
args = parentArgs.parse_args()

sm.setParams(0.12, 0.07, 0.90, 0.93, 0.20, 0.11)
sm.setSubPlots(2,1,"Pyrazine Cross Sections", 'Electron Energy (rydbergs)', 'Cross Section (bohr^2)')

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  totXss = kmat.getTotalDiscreteXS()
  totXss.setConversionEnergy(sm.eVs)
  totXss.writeAllToFile()
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()