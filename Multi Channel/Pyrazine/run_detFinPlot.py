from RatSMatWrap import *
import General.SimpPlot as sp
import argparse

parentArgs = argparse.ArgumentParser(description="Pyrazine Fit - Plot det(Fin)")
parentArgs.add_argument("plotStart_", help="Plot Start Energy", type=float)
parentArgs.add_argument("plotEnd_", help="Plot End Energy", type=float)
parentArgs.add_argument("plotComplex_", help="Plot Complex Energy Offset", type=float)
parentArgs.add_argument("steps_", help="Number of steps", type=int)
parentArgs.add_argument("subN_", help="Sub Set Number", type=int, nargs='?', default=None)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=None)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=None)
args = parentArgs.parse_args()

kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
vals = kmat.getFinDetRange(args.plotStart_, args.plotEnd_, args.plotComplex_, args.steps_)

sp.plotSingle("Pyrazine : det(Fin)", vals[0], [vals[1],vals[2]], 'Electron Energy (rydbergs)', 'det(Fin)', ['real','imag'])
