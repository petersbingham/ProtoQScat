import sys
sys.path.append("..")
from rfortmatreader import *

FILENAME = "fort.19"

NUMCHANNELS = 6
ENEFACTOR = sm.EFROMK_RYDBERGS
def getFileHandler(kCal, startIndex, endIndex):
    sysName = "Urazil_6chanEleastic_" + str(kCal) + "_" + str(startIndex)+ "_" + str(endIndex)
    return ResultFileHandler(sysName)

def readkMats(fileName):
    r = RfortMatReader(fileName, NUMCHANNELS)
    return r.readkMats()