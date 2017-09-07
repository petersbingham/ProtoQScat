PRECISION = 8
MIN = pow(10,-PRECISION)

def removeduplicateFloats(lst,comparator,accessor=None):
    def _getVal(lst,i,accessor):
        if accessor:
            return accessor(lst[i])
        else:
            return lst[i]
      
    if len(lst) > 0:
        iStart = 0
        while iStart<len(lst)-1:
            comVal = _getVal(lst,iStart,accessor)
            for i in range(iStart+1, len(lst)):
                if comparator.complexCompare(comVal,_getVal(lst,i,accessor)):
                    lst.pop(i)
                    break
                elif i == len(lst)-1:
                    iStart += 1          
          
class Compare:
    TYPE_ABSOLUTE = 0
    TYPE_PERCENT = 1
    def __init__(self, *args):
        if len(args) == 0:
            self.realType = self.TYPE_ABSOLUTE
            self.realVal = MIN
            self.imagType = self.TYPE_ABSOLUTE
            self.imagVal = MIN
        if len(args) == 1:
            self.realType = self.TYPE_ABSOLUTE
            self.realVal = args[0]
            self.imagType = self.TYPE_ABSOLUTE
            self.imagVal = args[0]
        elif len(args) == 2:
            self.realType = args[0]
            self.realVal = args[1]
            self.imagType = args[0]
            self.imagVal = args[1]
        elif len(args) == 4:
            self.realType = args[0]
            self.realVal = args[1]
            self.imagType = args[2]
            self.imagVal = args[3]
    
    def floatCompare(self, f1, f2, complex=False):
        if complex:
            return self._compare(f1, f2, self.imagType, self.imagVal)
        else:
            return self._compare(f1, f2, self.realType, self.realVal)
      
    def complexCompare(self, c1, c2):
        if not self.floatCompare(complex(c1).real, complex(c2).real, False):
            return False
        if not self.floatCompare(complex(c1).imag, complex(c2).imag, True):
            return False
        return True
    
    def _compare(self, f1, f2, type, val):
        if type == self.TYPE_ABSOLUTE:
            return self._absCmp(f1, f2, val)
        else:
            return self._relCmp(f1, f2, val)    
          
    def _absCmp(self, f1, f2, cmpVal):
        if abs(f1-f2) > cmpVal:
            return False
        return True
    
    def _relCmp(self, f1, f2, cmpVal):
        if abs(f1-f2)/f1*100 > cmpVal:
            return False
        return True


class RationalCompare:
    def __init__(self, zeroValue=0.0, distThreshold=None):
        self.distThreshold = distThreshold
        self.zeroValue = zeroValue

    def isClose(self, cval1, cval2, distThreshold=None):
        diff = self.getComplexDiff(cval1, cval2)
        return self.checkComplexDiff(diff, distThreshold)

    def getComplexDiff(self, cval1, cval2):
        realdiff = self._getDiff(cval2.real, cval1.real)
        imagdiff = self._getDiff(cval2.imag, cval1.imag)
        return realdiff + 1.0j*imagdiff

    def checkComplexDiff(self, cdiff, distThreshold=None):
        if distThreshold is None:
            distThreshold = self.distThreshold
        return cdiff.imag<distThreshold and cdiff.real<distThreshold

    def _getDiff(self, val1, val2):
        absVal1 = abs(val1)
        absVal2 = abs(val2)
        if absVal1<=self.zeroValue and absVal2<=self.zeroValue:
            return 0.0
        elif absVal1>self.zeroValue and absVal2>self.zeroValue:
            return self._calDiff(val1, val2)
        elif absVal1 > self.zeroValue:
            return self._calDiff(val1, self.zeroValue)
        else:
            return self._calDiff(self.zeroValue, val2)

    def _calDiff(self, val1, val2): #https://en.wikipedia.org/wiki/Relative_change_and_difference
        absVal1 = abs(val1)
        absVal2 = abs(val2)
        valMax = absVal2 if (absVal1 < absVal2) else absVal1
        return abs(val2-val1) / valMax


def absDiff(c1, c2):
    return abs(c1-c2)

def complexDiff(c1, c2):
    return c1-c2

def truncateFloat(dps, val):
    return float(("{0:."+str(dps)+"f}").format(val))

def truncateComplex(dps, val):
    return truncateFloat(dps,val.real) + truncateFloat(dps,val.imag)*1.0j

def getPermutations(possVals, numOf):
    if numOf == 1:
        perms = [[x] for x in possVals]
    else:
        perms = []
        nextPerms = getPermutations(possVals, numOf-1)
        for val in possVals:
            for nextPerm in nextPerms:
                perms.append([val]+nextPerm)
    return perms
   
def format_e(n):
    a = '%E' % n
    return a.split('E')[0].rstrip('0').rstrip('.') + 'e' + a.split('E')[1]         