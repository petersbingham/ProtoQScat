import numpy as np
import cmath
import matplotlib.pyplot as plt
import random

def decimate(mats, startIndex, endIndex, N):
    step = int((endIndex-startIndex) / (N-1))
    newMats = {}
    index = 0
    stepCnt = 0
    for ene in sorted(mats, key=lambda val: val.real):
        if index>endIndex:
            return None, None
        if index>=startIndex:
            if stepCnt == 0:
                newMats[ene] = mats[ene]
            stepCnt += 1
            if stepCnt == step:
                stepCnt = 0
        if len(newMats) == N:
            break
        index += 1
    return newMats, step, index

def getSfromKmatrices(kmats, numChannels):
    smats = {}
    for ene in sorted(kmats.keys()):
        mat = getSfromKmatrix(kmats, numChannels, ene)
        smats[ene] = mat
    return smats

def getSfromKmatrix(kmats, numChannels, ene):
    num = np.identity(numChannels) + 1.0j*kmats[ene]
    denum = np.identity(numChannels) - 1.0j*kmats[ene]
    S = np.dot(num, np.linalg.inv(denum))
    return S

REDUCED_MASS = 1.0     
K_POS = "kPos"
K_SIGN = "kSign"
K_ROT = "kRot"

class kCalculator:
    def __init__(self, thresholds, eneFactor, type=K_POS, channelSigns=None, invertChannel=False):  #invertChannel is for test purposes.
        self.thresholds = thresholds
        self.eneFactor = eneFactor
        self.type = type
        self.channelSigns = channelSigns
        self.invertChannel = invertChannel
    def __str__(self):
        if self.channelSigns is None:
            return self.type
        else:
            return str(self.channelSigns)
    def k(self, ch, ene):
        if self.type == K_POS:
            return self.kpos(ch, ene)
        elif self.type == K_SIGN:
            return self.ksign(ch, ene)
        elif self.type == K_ROT:
            return self.krot(ch, ene)
    def kpos(self, ch, ene):
        return cmath.sqrt(self._getValue(ch, ene))
    def ksign(self, ch, ene):
        if self.invertChannel:
            ch = len(self.thresholds)-1 - ch
        mult = 1.0
        if self.channelSigns is not None:
            mult = self.channelSigns[ch]
        return mult * self.kpos(ch, ene)
    def krot(self, ch, ene):
        k = self.kpos(ch, ene)
        absolute, argument = cmath.polar(k) 
        return absolute * cmath.exp(1j*(argument+self.getPhase(ch, ene)))
    def _getValue(self, ch, ene):
        return self.eneFactor*REDUCED_MASS*(ene - self.thresholds[ch])
    def getPhase(self, ch, ene):
        if ene.real <= self.thresholds[ch]:
            return 0.0
        else:
            return cmath.pi

class MatException(Exception):
    def __init__(self, string):
        self.string = string
    def __str__(self):
        return "Matrix Error: " + self.string
   
class matSequence:
    def __init__(self, title=None, colourCycle=None):
        if colourCycle is None:
            colourCycle = ['red', 'green', 'blue', 'purple']
        self.items = {}
        self.title = title
        self.marker = False
        self.legPrefix = ""
        self.colourCycle = colourCycle
    
    def __setitem__(self, key, item):
        self.size = item.size
        self.items[key] = item
    
    def writeAllToFile(self):
        f = open("res.txt","w")
        for x in sorted(self.items.keys(), key=lambda val: val.real):
            eneStr = "{:1.6f}".format(x.real) + "+{:1.6f}i".format(x.imag)
            f.write(eneStr + ":\n" + str(self.items[x])+"\n\n")
        f.close()
    
    def setDetails(self, legPrefix, colourCycle):
        self.legPrefix = legPrefix
        self.colourCycle = colourCycle
      
    def useMarker(self):
        self.marker = True
    
    def _init(self, ax):
        if self.title is not None:
            if ax==plt:
                fig = plt.figure()
                fig.suptitle(self.title)
            else:
                ax.set_title(self.title)
                ax.title.set_fontsize(10)
    
    def plot(self, m, n, logx=False, logy=False, imag=False, ax=plt):
        self._init(ax)
        l1,s1 = self._plot(m, n, logx, logy, imag, ax)
        if ax==plt and self.title is not None:
            ax.legend([l1], [s1])
            ax.draw()
        return (l1,s1)
    
    def _plot(self, m, n, logx, logy, imag, ax):
        items = self._convert()
        xs = np.ndarray((len(items),), dtype=float)
        ys = np.ndarray((len(items),), dtype=float)
        i = 0
        for x in sorted(items.keys(), key=lambda val: val.real):
            xs[i] = x.real
            if not imag:
                ys[i] = items[x][m,n].real
            else:
                ys[i] = items[x][m,n].imag
            i+=1
        if self.marker:
            ma = "x"
            li = "None"
        else:
            ma = None
            li = '-'
        if logx and logy:
            lne, = ax.loglog(xs, ys, linestyle=li, marker=ma, basex=10)
        elif logx:
            lne, = ax.semilogx(xs, ys, linestyle=li, marker=ma, basex=10)
        elif logy:
            lne, = ax.semilogy(xs, ys, linestyle=li, marker=ma, basey=10)
        else:
            lne, = ax.plot(xs, ys)
        return (lne, self.legPrefix+": "+str(m)+","+str(n))
        
    def plotRow(self, m, logx=False, logy=False, imag=False, ax=plt):
        self._init(ax)
        return self._plotMany(self._plotRow, logx, logy, imag, ax, m)
        
    def plotAll(self, logx=False, logy=False, imag=False, ax=plt):
        self._init(ax)
        return self._plotMany(self._plotAll, logx, logy, imag, ax)
      
    def _plotMany(self, plotFun, logx, logy, imag, ax, *args):
        plt.gca().set_color_cycle(self.colourCycle)
        lines, strings = plotFun(logx, logy, imag, ax, *args)
        if ax==plt and self.title is not None:
            plt.legend(lines, strings, prop={'size':9})
            plt.draw()
        return (lines, strings)
    
    def _plotAll(self, logx, logy, imag, ax):  
        size = self._getSize()
        lines = []
        strings = []
        for m in range(size):
            for n in range(size):
                l,s = self._plot(m, n, logx, logy, imag, ax)
                lines.append(l)
                strings.append(s)
        return (lines, strings)
    
    def _plotRow(self, logx, logy, imag, ax, *args):  
        size = self._getSize()
        lines = []
        strings = []
        for n in range(size):
            l,s = self._plot(args[0], n, logx, logy, imag, ax)
            lines.append(l)
            strings.append(s)
        return (lines, strings)
    
    def _getSize(self):
        key = random.choice(self.items.keys())
        return self.items[key].shape[0] 
    
    def _convert(self):
        return self.items

### Functionality for plotting several matrices.

subplots = None
subplotFig = None
numSubplotRows = 1
numSubplotCols = 1
subplotRowCnt = 1
subplotColCnt = 0

def setSubPlots(numSubplotRows_, numSubplotCols_, title=None, xlabel="", ylabel=""):
    global subplots
    global numSubplotRows
    global numSubplotCols
    global subplotFig
    numSubplotRows = numSubplotRows_
    numSubplotCols = numSubplotCols_
    subplotFig, subplots = plt.subplots(numSubplotRows, numSubplotCols, sharex='col', sharey='row', squeeze=False)
    if title is not None:
        plt.suptitle(title)
    subplotFig.text(0.5, 0.04, xlabel, ha='center', va='center')
    subplotFig.text(0.06, 0.5, ylabel, ha='center', va='center', rotation='vertical')

def _getAxis():
    global subplots
    if subplots is not None:
        global numSubplotRows
        global numSubplotCols
        global subplotRowCnt
        global subplotColCnt
        if subplotColCnt==numSubplotCols:
            subplotColCnt = 1
            subplotRowCnt += 1
        else:
            subplotColCnt += 1
        return subplots[subplotRowCnt-1][subplotColCnt-1]
    return plt

def plot(type, matSequences, m, n, imag=False, title=None, xlabel="", ylabel=""):
    _plot(type, _plot, matSequences, title, xlabel, ylabel, m, n, imag)  
    
def plotAll(type, matSequences, imag=False, title=None, xlabel="", ylabel=""):
    _plot(_plotAll, type, matSequences, title, xlabel, ylabel, imag)  
    
def plotSingle(type, matSequences, m, n, imag=False, title=None, xlabel="", ylabel=""):
    _plot(_plotSingle, type, matSequences, title, xlabel, ylabel, imag, m, n) 
    
def plotRow(type, matSequences, m, imag=False, title=None, xlabel="", ylabel=""):
    _plot(_plotRow, type, matSequences, title, xlabel, ylabel, imag, m) 

def _plot(plotFun, type, matSequences, title, xlabel, ylabel, imag, *args):
    global subplots
    global totNumSubplots
    global numSubplots
    if title is not None:
        if subplots is None:
            fig = plt.figure()
            fig.suptitle(title)
        else:
            for matSequence in matSequences:
                matSequence.title = title
    lines, strings = plotFun(type, matSequences, imag, *args)
    if subplots is None:
        plt.legend(lines, strings, prop={'size':9})
        if xlabel != "":
            plt.xlabel(xlabel, fontsize=12)
        if ylabel != "":
            plt.ylabel(ylabel, fontsize=12)
    else:
        plt.figlegend(lines, strings, prop={'size':9}, loc = 'upper left', ncol=1, labelspacing=0. )
    if subplots is None or (numSubplotRows==subplotRowCnt and numSubplotCols==subplotColCnt):
        plt.draw()
        plt.show()

def _plotAll(type, matSequences, imag):
    lines = []
    strings = []
    for matSequence in matSequences:
        ax = _getAxis()
        if type == "Lin":
            l,s = matSequence.plotAll(False, False, imag, ax)
        elif type=="Log":
            l,s = matSequence.plotAll(False, True, imag, ax)
        elif type=="LogLog":
            l,s = matSequence.plotAll(True, True, imag, ax)
        lines += l
        strings += s
    return (lines, strings)

def _plotSingle(type, matSequences, imag, *args):
    lines = []
    strings = []
    for matSequence in matSequences:
        ax = _getAxis()
        if type == "Lin":
            l,s = matSequence.plot(args[0], args[1], False, False, imag, ax)
        elif type=="Log":
            l,s = matSequence.plot(args[0], args[1], False, True, imag, ax)
        elif type=="LogLog":
            l,s = matSequence.plot(args[0], args[1], True, True, imag, ax)
        lines.append(l)
        strings.append(s)
    return (lines, strings)

def _plotRow(type, matSequences, imag, *args):
    lines = []
    strings = []
    for matSequence in matSequences:
        ax = _getAxis()
        if type == "Lin":
            l,s = matSequence.plotRow(args[0], False, False, imag, ax)
        elif type=="Log":
            l,s = matSequence.plotRow(args[0], False, True, imag, ax)
        elif type=="LogLog":
            l,s = matSequence.plotRow(args[0], True, True, imag, ax)
        lines += l
        strings += s
    return (lines, strings)
  
class mat:
    PADDING = 3
    def __init__(self, size, precision):
        self.size = size
        self.precision = precision
        self.min = pow(10,-self.precision)
        self.formatstr = "{:+."+str(self.precision)+"f}"
        
    def __getitem__(self, i):
        return self._getRow(i)
      
    def getMatrix(self):
        mlist = []
        for m in range(self.size):
            rlist = []
            for n in range(self.size):
                rlist.append(complex(self[m][n]))
            mlist.append(rlist)
        return np.matrix(mlist)
    
    def __str__(self):
        isImag = self._isImag()
        
        maxLen = 0
        for m in range(self.size-1):
            for n in range(self.size-1):
                eLen = len(self._getFormattedStr(self[m][n], isImag))
                if eLen > maxLen:
                    maxLen = eLen 
        s = ""     
        for m in range(self.size):
            for n in range(self.size):
                s += self._padStr(self[m][n], isImag, maxLen)
            s += "\n"   
        return s
  
    def _padStr(self, value, isImag, maxLen=None):
        theStr = self._getFormattedStr(value, isImag)
        if maxLen is not None:
            return theStr.ljust(maxLen + self.PADDING)  
        else: 
            return theStr
      
    def _isImag(self):
        for m in range(0,self.size):
            for n in range(0,self.size):
                if abs(float(complex(self[m][n]).imag)) > self.min:
                    return True
        return False
    
    def _getFormattedStr(self, value, isImag):
        if isImag:
            return "("+self.formatstr.format(complex(value).real) + self.formatstr.format(complex(value).imag)+"i)"
        else:
            return self.formatstr.format(complex(value).real)
    