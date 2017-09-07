import sys
sys.path.append("../utilities")
import scattering.matrices as sm

class Decimator():
    def __init__(self, startIndex, endIndex, offset, resultFileHandler):
        self.startIndex = startIndex
        self.endIndex = endIndex
        self.offset = offset
        self.resultFileHandler = resultFileHandler
                 
    def decimate(self, sMats, N):
        actualStartIndex = self.startIndex+self.offset
        sMats, step, actualEndIndex, startEne, endEne = sm.decimate(sMats, actualStartIndex, self.endIndex+self.offset, N)
        self.resultFileHandler.setDeciminationInfo(actualStartIndex, actualEndIndex)
        try:
            decStr = "N=%d, Emin=%d(%f), Emax=%d(%f), step=%d" % (N,actualStartIndex,startEne,actualEndIndex,endEne,step)
        except TypeError:
            decStr = "N=%d, Emin=%d(%f+%fi), Emax=%d(%f+%fi), step=%d" % (N,actualStartIndex,startEne.real,startEne.imag,actualEndIndex,endEne.real,endEne.imag,step)
        print "Decimation:"
        print "  "+decStr
        return sMats, decStr  