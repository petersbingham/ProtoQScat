import sys
sys.path.append("..")
from rfortmatreader import *

FILENAME = "fort.19"

NUMCHANNELS = 3
ENEFACTOR = sm.EFROMK_RYDBERGS
def getFileHandler(kCal, startIndex, endIndex):
    sysName = "Pyrazine" + "_" + str(kCal) + "_" + str(startIndex)+ "_" + str(endIndex)
    return ResultFileHandler(sysName)

def readkMats(fileName):
    r = RfortMatReader(fileName, NUMCHANNELS)
    return r.readkMats()