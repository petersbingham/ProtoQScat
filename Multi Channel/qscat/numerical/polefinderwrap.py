import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')
sys.path.insert(0,os.getcwd()) #We assume that the specific kreader and description (below) will be here.

from elastic3chanreader import *
from description import *
from polefinder import *

MODE_ALLSIGNS_DOUBLE = 0
MODE_ALLSIGNS_INC = 1
MODE_ROT_DOUBLE = 2
MODE_ROT_INC = 3
MODE_POS_DOUBLE = 4
MODE_POS_INC = 5

import argparse
parentArgs = argparse.ArgumentParser(description="Pyrazine Fit - Pole find")
parentArgs.add_argument("startIndex_", help="Start Index", type=int)
parentArgs.add_argument("endIndex_", help="End Index", type=int)
parentArgs.add_argument("offset_", help="Offset", type=int)
parentArgs.add_argument("mode_", help="Mode", type=int)
parentArgs.add_argument("cfSteps_", help="Compare Steps", type=int)
parentArgs.add_argument("distFactor_", help="Distinguish Factor", type=float)
parentArgs.add_argument("zeroValExp_", help="Zero Value Precision", type=int)
parentArgs.add_argument("cmpValue_", help="Compare Value", type=complex, nargs='?', default=None)
args = parentArgs.parse_args()

def _doPoleFind(kCal, mode, dirName):
    kmats = readkMats("./fort.19")
    resultFileHandler = getFileHandler(kCal, args.startIndex_, args.endIndex_)
    PoleFinder(sm.getSfromKmatrices(kmats,NUMCHANNELS), kCal, resultFileHandler, ENEFACTOR, args.startIndex_, args.endIndex_, args.offset_, args.distFactor_, args.cfSteps_, args.cmpValue_, mode)

def _polesForAllSigns(mode, dirName):
    for i in [1.0,-1.0]:
      for j in [1.0,-1.0]:
        for k in [1.0,-1.0]:
          ksigns = [i,j,k]
          kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=ksigns, eneFactor=ENEFACTOR)
          _doPoleFind(kCal, mode, dirName)

def _polesForRot(mode, dirName):
    kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_ROT, eneFactor=ENEFACTOR)
    _doPoleFind(kCal, mode, dirName)

def _polesForPos(mode, dirName):
    kCal = sm.kCalculator(THRESHOLDS, LS, ktype=sm.K_SIGN, ksigns=[1.0,1.0,1.0], eneFactor=ENEFACTOR)
    _doPoleFind(kCal, mode, dirName)

if args.mode_ == MODE_ALLSIGNS_DOUBLE:
    _polesForAllSigns(DOUBLE_N, "AllSigns Double")
elif args.mode_ == MODE_ALLSIGNS_INC:
    _polesForAllSigns(INC_N, "AllSigns Inc")
elif args.mode_ == MODE_ROT_DOUBLE:
    _polesForRot(DOUBLE_N, "Rot Double")
elif args.mode_ == MODE_ROT_INC:
    _polesForRot(INC_N, "Rot Inc")
elif args.mode_ == MODE_POS_DOUBLE:
    _polesForPos(DOUBLE_N, "Pos Double")
elif args.mode_ == MODE_POS_INC:
    _polesForPos(INC_N, "Pos Inc")