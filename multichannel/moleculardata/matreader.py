import sys
sys.path.append("..")
from rfortmatreader import *
from sysdesc import *

MASSMULT = sm.MASSMULT_RYDBERGS
def getFileHandler(kCal, startIndex, endIndex):
    sysName = ARCHIVE_BASE_STR + "_" + str(kCal) + "_" + str(startIndex)+ "_" + str(endIndex)
    return ResultFileHandler(sysName)

def readkMats(fileName):
    r = RfortMatReader(fileName, NUMCHANNELS)
    return r.readkMats()