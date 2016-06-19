import matplotlib
import matplotlib.pyplot as plt
import numpy as np
import math
import csv

LEFT = 0.12
BOTTOM = 0.1
RIGHT = 0.9
TOP = 0.9
WSPACE = 0.2
HSPACE = 0.2

FIGSZ_W = 8
FIGSZ_H = 6
  
class StaticPlot:
    DPI = 80
  
    def __init__(self, title, rows=1, cols=1, drawAxes=False):
        self.fig = plt.figure(figsize=(FIGSZ_W,FIGSZ_H), dpi=self.DPI)
        self.fig.suptitle(title)
        self.rows = rows
        self.cols = cols
        self.drawAxes = drawAxes
        self.lines = []
        self.legends = []
        self.axisConfig = []
        self.fig.subplots_adjust(left=LEFT, bottom=BOTTOM, right=RIGHT, top=TOP, wspace=WSPACE, hspace=HSPACE)
  
    def addPlot(self, xlabel, ylabel, logx=False, logy=False):
        plotNum = len(self.axisConfig)+1
        self.fig.add_subplot(self.rows,self.cols,plotNum)
        plt.xlabel(xlabel, fontsize=12)
        plt.ylabel(ylabel, fontsize=12)
        self.axisConfig.append([logx,logy])
    
    def addLine(self, plotNum, xs, ys, legend=None, markerSz=None, markWithLine=False):
        self._addData(plotNum, xs, ys, legend, markerSz, False, markWithLine)
    
    def addScat(self, plotNum, xs, ys, legend=None, markerSz=None):
        self._addData(plotNum, xs, ys, legend, markerSz, True, False)
  
    def _addData(self, plotNum, xs, ys, legend, markerSz, scatter, markWithLine):
        xs, ys = self._convertValues(xs, ys)
        if scatter:
            l = self._addScatType(plotNum, xs, ys, markerSz) 
        else:
            l = self._addLineType(plotNum, xs, ys, markerSz, markWithLine)  
        self.lines.append(l)
        if legend is not None:
            self.legends.append(legend)
  
    def _convertValues(self, xs, ys):
        if type(xs) is list:
            xs = np.ndarray((len(xs),), buffer=np.array(xs)) #Contents of array need to be floats or can get: "TypeError: buffer is too small for requested array"
        if type(ys) is list:
            ys = np.ndarray((len(ys),), buffer=np.array(ys))
        return (xs,ys)
  
    def _addScatType(self, plotNum, xs, ys, markerSz):
        plt.scatter(xs,ys,color='blue',s=5,edgecolor='none')
  
    def _addLineType(self, plotNum, xs, ys, markerSz, markWithLine):
        kwargs = {'basex':10}
        if markerSz:
            kwargs = {'basex':10, 'linestyle':'None' if not markWithLine else 'solid', 'marker':'o', 'markerfacecolor':'r', 'markersize':markerSz}
        if self.axisConfig[plotNum-1][0] and self.axisConfig[plotNum-1][1]:
            l, = plt.loglog(xs, ys, **kwargs)
        elif self.axisConfig[plotNum-1][0]:
            l, = plt.semilogx(xs, ys, **kwargs)
        elif self.axisConfig[plotNum-1][1]:
            kwargs.pop('basex')
            l, = plt.semilogy(xs, ys, **kwargs)
        else:
            kwargs.pop('basex')
            #print str(len(xs)) + " " + str(len(ys))
            l, = plt.plot(xs, ys, **kwargs)
        return l
  
    def reveal(self):  
        if self.drawAxes:
            plt.axhline(0, color='black')
            plt.axvline(0, color='black')
        if len(self.legends) > 0:
            plt.legend(self.lines, self.legends, prop={'size':9})
        axes = plt.gca()
        if xlim is not None:
            axes.set_xlim(xlim)
        if ylim is not None:
            axes.set_ylim(ylim)
        plt.draw()
        plt.show()
    
    def save(self, path):
        self.fig.savefig(path, dpi=self.DPI)
  
def setSubPlotParameters(left=None, bottom=None, right=None, top=None, wspace=None, hspace=None):
    global LEFT
    global BOTTOM
    global RIGHT
    global TOP
    global WSPACE
    global HSPACE
    if left: LEFT = left
    if bottom: BOTTOM = bottom
    if right: RIGHT = right
    if top: TOP = top
    if wspace: WSPACE = wspace
    if hspace: HSPACE = hspace

def setImgSize(width, height):
    global FIGSZ_W
    global FIGSZ_H
    FIGSZ_W = width
    FIGSZ_H = height

xlim = None
ylim = None
def setExtents(xlim_, ylim_):
    global xlim
    global ylim
    xlim = xlim_
    ylim = ylim_
    
def turnOffColourCycle():
    matplotlib.rcParams['axes.color_cycle'] = ['black']

def _getDataFromCSV(csvpath):
    xs = []
    yss = None
    with open(csvpath, 'rb') as csvfile:
        csvReader = csv.reader(csvfile, delimiter=',')
        for row in csvReader:
            xs.append(float(row[0]))
            if yss is None:
                yss = [[] for i in range(0,len(row[1:]))]
            for i in range(0,len(row[1:])):
                yss[i].append(float(row[1+i]))
    return xs, yss
            
def plotSingleFromCSV(csvpath,title,xlabel="",ylabel="",legends=None,logx=False,logy=False,markerSz=None,path=None,markWithLine=False,drawAxes=False):
    xs,yss = _getDataFromCSV(csvpath)
    plotSingle(title,xs,yss,xlabel,ylabel,legends,logx,logy,markerSz,markWithLine,path,drawAxes)
    
def plotSingle(title,xs,yss,xlabel="",ylabel="",legends=None,logx=False,logy=False,markerSz=None,path=None,markWithLine=False,drawAxes=False):
    _plot(title,xs,yss,xlabel,ylabel,legends,logx,logy,markerSz,markWithLine,path,drawAxes)
  
def plotSingle2(title,xss,yss,xlabel="",ylabel="",legends=None,logx=False,logy=False,markerSz=None,path=None,markWithLine=False,drawAxes=False):
    _plot2(title,xss,yss,xlabel,ylabel,legends,logx,logy,markerSz,markWithLine,path,drawAxes)

def _plot(title,xs,yss,xlabel,ylabel,legends,logx,logy,markerSz,markWithLine,path,drawAxes):
    p = StaticPlot(title,drawAxes=drawAxes)
    p.addPlot(xlabel, ylabel, logx, logy)
    for i in range(len(yss)):
        if legends is not None:
            p.addLine(1,xs,yss[i],legends[i],markerSz,markWithLine)
        else:
            p.addLine(1,xs,yss[i],None,markerSz,markWithLine)
    p.reveal()
    if path is not None:
        p.save(path)

def _plot2(title,xss,yss,xlabel,ylabel,legends,logx,logy,markerSz,markWithLine,path,drawAxes):
    p = StaticPlot(title,drawAxes=drawAxes)
    p.addPlot(xlabel, ylabel, logx, logy)
    for i in range(len(xss)):
        if legends is not None:
            p.addLine(1,xss[i],yss[i],legends[i],markerSz,markWithLine)
        else:
            p.addLine(1,xss[i],yss[i],None,markerSz,markWithLine)
    p.reveal()
    if path is not None:
        p.save(path)
  
def plotScat(title,xs,ys,xlabel="",ylabel="",logx=False,logy=False,legend=None,markerSz=None,path=None,drawAxes=False):
    p = StaticPlot(title,drawAxes=drawAxes)
    p.addPlot(xlabel, ylabel, logx, logy)
    p.addScat(1,xs,ys,legend,markerSz)
    p.reveal()
    if path is not None:
        p.save(path)  
  
if __name__ == '__main__':
    numValues = 400
    def _getWave(xs, ys, phase, amp):  
        x = 0.0
        dx = 4.0*np.pi / float(numValues)
        for i in range(numValues):
            xs[i] = x
            ys[i] = amp*math.sin(x+phase)
            x += dx
        return (xs,ys)  
    
    def _getWaveNumpy(phase, amp):
        xs = np.ndarray((numValues,), dtype=float)
        ys = np.ndarray((numValues,), dtype=float)
        return _getWave(xs, ys, phase, amp)
    
    def _getWaveList(phase, amp):
        xs = [None]*numValues
        ys = [None]*numValues
        return _getWave(xs, ys, phase, amp)
  
    def _getLog():
        xs = np.ndarray((numValues,), dtype=float)
        ys = np.ndarray((numValues,), dtype=float)
        x = 0.0
        dx = 10.0 / float(numValues)
        for i in range(numValues):
            xs[i] = x
            ys[i] = math.pow(2,x)
            x += dx
        return (xs,ys)  
    
    sp = StaticPlot("Test",1,2)
      
    sp.addPlot("T1","T2")
    xs,ys = _getWaveNumpy(0,1)
    sp.addLine(1,xs,ys,"L1",True)
    
    xs,ys = _getWaveList(np.pi/2,2)
    sp.addLine(1,xs,ys,"L2")
    
    sp.addPlot("T1","T2", logy=True)
    xs,ys = _getLog()
    sp.addLine(2,xs,ys,"L2")
    
    sp.reveal()

