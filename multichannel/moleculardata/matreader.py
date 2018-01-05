import sys
sys.path.append("..")
from rfortmatreader import *
from sysdesc import *

MASSMULT = sm.MASSMULT_RYDBERGS
def getFileHandler(kCal, startIndex, endIndex, fromEnd=False):
    sysName = ARCHIVE_BASE_STR + "_" + str(kCal) + "_" + str(startIndex)+ "_" + str(endIndex)
    if fromEnd:
        sysName += "_EndLock"
    return ResultFileHandler(sysName)

def readkMats(fileName):
    r = RfortMatReader(fileName)
    return r.readkMats()