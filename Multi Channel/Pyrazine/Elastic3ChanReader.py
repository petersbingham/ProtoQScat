import sys
sys.path.append("..")
sys.path.append("../..")
sys.path.append("../../Utilities")
sys.path.append("../../../Utilities")
import Scattering.Matrices as sm
from General.QSType import *

NUMCHANNELS = 3
ENEFACTOR = sm.EFROMK_RYDBERGS

def getFitName(kCal, startIndex, endIndex):
    return "Pyrazine" + "_" + str(kCal) + "_" + str(startIndex)+ "_" + str(endIndex)

def readkMats(fileName):
    numChannels = 3
    kmats = {}
    with open(fileName, 'rb') as file:
      linNum = 0
      for i in range(0,5):
        linNum += 1
        file.readline()
      ene = None
      for line in file:
        linNum += 1
        nums = line.split()
        if len(nums) == 4 and nums[0]=="3":
          ene = _num(nums[3])
          kmats[ene] = QSsqZeros(numChannels)
        else:
          if len(line) == 81:
            kmats[ene][0,0] = _num(line[0:20])
            kmats[ene][1,0] = _num(line[21:40])
            kmats[ene][1,1] = _num(line[41:60])
            kmats[ene][2,0] = _num(line[61:80])
            kmats[ene][0,1] = kmats[ene][1,0]  
            kmats[ene][0,2] = kmats[ene][2,0]
          elif len(line) == 41:
            kmats[ene][2,1] = _num(line[0:20])
            kmats[ene][2,2] = _num(line[21:40])
            kmats[ene][1,2] = kmats[ene][2,1]
          elif len(line) != 0:
            raise Exception("Line " + str(linNum) + " bad: " + str(line) + "  Len: " + str(len(line)))
    return kmats

def _num(string):
  return float(string.replace("D", "E").replace(" ", ""))