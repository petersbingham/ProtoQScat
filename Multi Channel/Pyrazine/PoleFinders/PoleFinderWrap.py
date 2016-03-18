import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/.')
sys.path.insert(0,base+'/..')
sys.path.insert(0,base+'/../..')
sys.path.insert(0,base+'/../../../Utilities')
from PoleFinder import *
from Elastic3ChanReader import *


import argparse
parentArgs = argparse.ArgumentParser(description="Pyrazine Fit - Pole find")
parentArgs.add_argument("startIndex_", help="Start Index", type=int)
parentArgs.add_argument("endIndex_", help="End Index", type=int)
parentArgs.add_argument("offset_", help="Offset", type=int)
parentArgs.add_argument("distFactor_", help="Distinguish Factor", type=float)
parentArgs.add_argument("cmpValue_", help="Compare Value", type=complex, nargs='?', default=None)
args = parentArgs.parse_args()

for i in [1.0,-1.0]:
  for j in [1.0,-1.0]:
    for k in [1.0,-1.0]:
      ksigns = [i,j,k]
      kCal = sm.kCalculator([0.0,0.0,0.0], EFROMK_CONVERSIONFACTOR, sm.K_SIGN, ksigns)
      kmats = readkMats("../fort.19")
      fitName = getFitName(kCal, args.startIndex_, args.endIndex_)
      PoleFinder(sm.getSfromKmatrices(kmats,NUMCHANNELS), kCal, "./Results", fitName, EFROMK_CONVERSIONFACTOR, args.startIndex_, args.endIndex_, args.offset_, args.distFactor_, args.cmpValue_)