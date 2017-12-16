import sys
import os
import __main__

sys.path.append(".")
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../utilities')
sys.path.insert(0,base+'../..')

if __name__ != '__main__': #Can test without these imports
    import scattering.matrices as sm
    import general.type_wrap as tw
    from ratsmat.resultfilehandler import *

NUMFIELDCHARS=20

class RfortMatReader:
    def __init__(self, fileName):
        self.fileName = fileName

    def readkMats(self):
        self.kmats = {}
        with open(self.fileName, 'rb') as file:
            linNum = 0
            for i in range(0,5):
                linNum += 1
                file.readline()
            ene = None
            numCompleteLinesPerMat = 0
            numRemElements = 0
            for line in file:
                linNum += 1
                nums = filter(lambda x: x!="", line.split())
                if nums[0].find(".") == -1:
                    numChannels = int(nums[0])
                    numUniqueElements = self._aSum(numChannels)
                    numCompleteLinesPerMat = numUniqueElements / 4
                    numRemElements = numUniqueElements % 4
                    lineI = 0
                    self.cElement = 1
                    ene = self._num(nums[3])
                    self.kmats[ene] = tw.sqZeros(numChannels)
                else:
                    if lineI < numCompleteLinesPerMat:
                        self._checkLineLength(line, linNum, 4)
                        self._readLines(ene, line, 4)
                    elif numRemElements > 0:
                        self._checkLineLength(line, linNum, numRemElements)
                        self._readLines(ene, line, numRemElements)
                    lineI += 1
        self._flipCopyDiag()        
        return self.kmats
            
    def _checkLineLength(self, line, linNum, numElements):
        if len(line)!=numElements*NUMFIELDCHARS+1 and len(line)!=numElements*NUMFIELDCHARS+2: #Depending on new line sequence
            raise Exception("Line " + str(linNum) + " bad: " + str(line) + "  Len: " + str(len(line)) + "  Elements: " + str(numElements))
    
    def _readLines(self, ene, line, numElements):
        for i in range(numElements):
            indices = self._getIndices()
            self.kmats[ene][indices[0],indices[1]] = self._num(line[i*20:(i+1)*20])
                                
    def _getIndices(self):
        rows = 1
        cnt = 0
        while cnt < self.cElement:
            cnt += rows
            rows += 1
        ri = rows-2
        pSum = self._aSum(rows-2)
        ci = self.cElement - pSum - 1
        self.cElement += 1
        return (ri, ci)
        
    def _aSum(self, n):
        return n * (n+1) / 2
    
    def _num(self, string):
        return float(string.replace("D", "E").replace(" ", ""))
    
    def _flipCopyDiag(self):
        for ene in self.kmats:
            kmat = self.kmats[ene]
            numChannels = tw.shape(kmat)[0]
            numUniqueElements = self._aSum(numChannels) 
            self.cElement = 1
            for i in range(numUniqueElements):
                indices = self._getIndices()
                if indices[0] != indices[1]:
                    kmat[indices[1],indices[0]] = kmat[indices[0],indices[1]]
    
if __name__ == '__main__':
    r = RfortMatReader(None)
    for i in range(3):
        print r._getIndices()