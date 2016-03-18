import sys
import os
import copy
import numpy as np
from tabulate import tabulate
from sets import Set
from fileinput import filename
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../Utilities')
sys.path.insert(0,base+'/../..')
import General.Numerical as num
from GenDat import *
  
LVL_A = 0
LVL_V = 1
LVL_L = 2
LVL_KMIN = 3
LVL_KMAX = 4
LVL_DF = 5
LVL_N = 6

LVL_LAST = 6

def makePathStr(a,V,l,kmin,kmax,df,N,lockLvl=-1):
  if lockLvl==-1:
    return str(a) + "_"  + str(V) + "_" + str(l) + "_" + str(kmin) + "_" + str(kmax) + "_" + str(df) + "_" + str(N)
  elif lockLvl==LVL_A:
    return "_x_" + str(V) + "_" + str(l) + "_" + str(kmin) + "_" + str(kmax) + "_" + str(df) + "_" + str(N)
  elif lockLvl==LVL_V:
    return str(a) + "_x_" + str(l) + "_" + str(kmin) + "_" + str(kmax) + "_" + str(df) + "_" + str(N)
  elif lockLvl==LVL_L:
    return str(a) + "_"  + str(V) + "_x_" + str(kmin) + "_" + str(kmax) + "_" + str(df) + "_" + str(N)
  elif lockLvl==LVL_KMIN:
    return str(a) + "_"  + str(V) + "_" + str(l) + "_x_" + str(kmax) + "_" + str(df) + "_" + str(N)
  elif lockLvl==LVL_KMAX:
    return str(a) + "_"  + str(V) + "_" + str(l) + "_" + str(kmin) + "_x_" + str(df) + "_" + str(N)
  elif lockLvl==LVL_DF:
    return str(a) + "_"  + str(V) + "_" + str(l) + "_" + str(kmin) + "_" + str(kmax) + "_x_" + str(N)
  elif lockLvl==LVL_N:
    return str(a) + "_"  + str(V) + "_" + str(l) + "_" + str(kmin) + "_" + str(kmax) + "_" + str(df) + "_x_"

def makePathStrFromLst(lst,lockLvl=-1):
  return makePathStr(lst[0],lst[1],lst[2],lst[3],lst[4],lst[5],lst[6],lockLvl)

################### Init Functions ###################
    
class Dat:
  def __init__(self):
    self.data = {}
    self.comparator = num.Compare(0.001)
    self.p = Printer()
    
  def setComparator(self, comparator):
    self.comparator = comparator

################### Data Build Functions ###################
    
  def addRoot(self,path,k,E,prevRootVal=None,diff=None):
    end = self._addPathToSelfData(path)
    end[INDEX_ROOT].append(Root(k,E,prevRootVal,diff))
  
  def addPole(self,path,k,E,isNew,isLost):
    end = self._addPathToSelfData(path)
    end[INDEX_POLE].append(Pole(k,E,isNew,isLost))
  
  def _addPathToSelfData(self,path):
    return self._addPath(self.data,path[LVL_A],path[LVL_V],path[LVL_L],path[LVL_KMIN],path[LVL_KMAX],path[LVL_DF],path[LVL_N])
  
  def _addPath(self,data,a,V,l,kmin,kmax,df,N):
    if a not in data:
      data[a] = {}
    if V not in data[a]:
      data[a][V] = {}
    if l not in data[a][V]:
      data[a][V][l] = {}
    if kmin not in data[a][V][l]:
      data[a][V][l][kmin] = {}
    if kmax not in data[a][V][l][kmin]:
      data[a][V][l][kmin][kmax] = {}
    if df not in data[a][V][l][kmin][kmax]:
      data[a][V][l][kmin][kmax][df] = {}
    if N not in data[a][V][l][kmin][kmax][df]:
      data[a][V][l][kmin][kmax][df][N] = [[],[],[]]
    end = data[a][V][l][kmin][kmax][df][N]
    return end
  
####################### Accessor #######################
  
  def getData(self):
    return self.data
  
######################### Info #########################
    
  def getRange(self, searchLvl):
    range = [None,None]
    self._R_getRange(self.data, searchLvl, range, 0)
    return range
   
  def _R_getRange(self, data, searchLvl, range, lvlCnt):
      for val in data:
        if lvlCnt==searchLvl:
          if range[0] is None or val < range[0]:
            range[0] = val
          if range[1] is None or val > range[1]:
            range[1] = val
        else:
          self._R_getRange(data[val], searchLvl, range, lvlCnt+1)
    
################### Search Functions ###################
    
  def searchkForRoot(self, searchRoot, lvl=None):
    return self._searchPath(searchRoot, lvl, INDEX_ROOT, self._getk)
  def searchEForRoot(self, searchRoot, lvl=None):
    return self._searchPath(searchRoot, lvl, INDEX_ROOT, self._getE)
  def searchkForPole(self, searchPole, lvl=None):
    return self._searchPath(searchPole, lvl, INDEX_POLE, self._getk)
  def searchEForPole(self, searchPole, lvl=None):
    return self._searchPath(searchPole, lvl, INDEX_POLE, self._getE)

  def listSearch(self, searchFun, searchVals):
    results = {}
    for val in searchVals:
        results[val] = searchFun(val)
    return results

  def dictSearch(self, searchFun, searchVals, lvl):
    results = {}
    for val in searchVals.items():
        results[val[1]] = searchFun(val, lvl)
    return results
    
  def _searchPath(self, searchValue, lvl, index, fun):
    data = {}
    cmpLvl = None
    for a in self.data:
      if lvl is not None and lvl==LVL_A: cmpLvl=a
      for V in self.data[a]:
        if lvl is not None and lvl==LVL_V: cmpLvl=V
        for l in self.data[a][V]:
          if lvl is not None and lvl==LVL_L: cmpLvl=l
          for kmin in self.data[a][V][l]:
            if lvl is not None and lvl==LVL_KMIN: cmpLvl=kmin
            for kmax in self.data[a][V][l][kmin]:
              if lvl is not None and lvl==LVL_KMAX: cmpLvl=kmax
              for df in self.data[a][V][l][kmin][kmax]:
                if lvl is not None and lvl==LVL_DF: cmpLvl=df
                for N in self.data[a][V][l][kmin][kmax][df]:
                  if lvl is not None and lvl==LVL_N: cmpLvl=N
                  for value in self.data[a][V][l][kmin][kmax][df][N][index]:
                    cmpValue = fun(value)
                    fnd=False
                    if lvl is None:
                      if self._compareComplex(searchValue, cmpValue):
                        fnd = True
                    else:
                      if self._compareComplex(searchValue[0], cmpLvl) and self._compareComplex(searchValue[1], cmpValue):
                        fnd = True
                    if fnd:
                      end = self._addPath(data,a,V,l,kmin,kmax,df,N)
                      end[index].append(value)
    return data
  
  def _compareComplex(self, val1, val2):
    return self.comparator.complexCompare(val1, val2)

  def _getk(self, value):
    return value.k
  def _getE(self, value):
    return value.E

class Parser:                  
  def isN(self, line):
    return line[0] == 'N'

  def getN(self, line):
    return int(line[2:])

  def isRoot(self, line):
    return "Root_" in line
    
  def setRoot(self, d, line, path, N):
    k,E = self._extractNums(line)
    prevRootVal, diff = self._extractRootInfo(line)
    d.addRoot(path+[N],k,E,prevRootVal,diff)

  def isPole(self, line):
    return "Pole_" in line and "MiscPole_" not in line

  def setPole(self, d, line, path, N):
    new = POLE_STATUS_NEW in line
    lost = POLE_STATUS_LOST in line
    k,E = self._extractNums(line)
    d.addPole(path+[N],k,E,new,lost)
    
  def _extractNums(self, line):
    repLine = line.replace('i','j')
    numStrings = repLine.split("\t")
    kString = numStrings[0].split("=")[1]
    eString = numStrings[1].split("=")[1]
    return [complex(kString), complex(eString)]
  
  def _extractRootInfo(self, line):
    a = line.find("diff")
    if a!=-1:
      b = line.index("]",a)
      c = line.index("=",a)
      nums = line[a+5:b].split(" ")
      diff = complex(line[c+2:].replace('i','j'))
      return (int(nums[1]), diff)
    return (None,None)
  
class Filter:
  def __init__(self):
    self.comparator = num.Compare(0.0000001)
    
  def selectWhereIn(self,data,a_s=None,V_s=None,l_s=None,kmin_s=None,kmax_s=None,df_s=None,N_s=None):
    itData = copy.deepcopy(data)
    self._R_select(data, itData, 0, [a_s,V_s,l_s,kmin_s,kmax_s,df_s,N_s], lambda crit,val : val in crit)
  
  def selectWhereLessThan(self,data,a=None,V=None,l=None,kmin=None,kmax=None,df=None,N=None):
    itData = copy.deepcopy(data)
    self._R_select(data, itData, 0, [a,V,l,kmin,kmax,df,N], lambda crit,val : val < crit)
  
  def selectWhereGreaterThan(self,data,a=None,V=None,l=None,kmin=None,kmax=None,df=None,N=None):
    itData = copy.deepcopy(data)
    self._R_select(data, itData, 0, [a,V,l,kmin,kmax,df,N], lambda crit,val : val > crit)

  def _R_select(self, data, itData, lvlCnt, critList, fun):
    for val in itData:
      if critList[lvlCnt] is not None and not fun(critList[lvlCnt], val):
        data.pop(val)
      elif lvlCnt<LVL_LAST:
        self._R_select(data[val], itData[val], lvlCnt+1, critList, fun)
  
  def clean(self,data):
    self._removeDuplicates(data)

  def _removeDuplicates(self,data):
    for a in data:
      for V in data[a]:
        for l in data[a][V]:
          for kmin in data[a][V][l]:
            for kmax in data[a][V][l][kmin]:
              for df in data[a][V][l][kmin][kmax]:
                for N in data[a][V][l][kmin][kmax][df]:
                  for type in data[a][V][l][kmin][kmax][df][N]:
                    num.removeduplicateFloats(type,self.comparator,lambda x: x.k)              

class Printer:
  def __init__(self, printReal=True, printImag=True):
    self.setPrintConfig(printReal, printImag)

  def setPrintConfig(self, printReal, printImag):
    self.printReal = printReal
    self.printImag = printImag
    
  def listParameters(self, data, ignoreN=True):
    table = [["a","V","l","kmin","kmax","df"]]
    if not ignoreN:
      table[0].append("N")
    row = []
    for a in sorted(data):
      self._addRow(row, a, LVL_A)
      for V in sorted(data[a]):
        self._addRow(row, V, LVL_V)
        for l in sorted(data[a][V]):
          self._addRow(row, l, LVL_L)
          for kmin in sorted(data[a][V][l]):
            self._addRow(row, kmin, LVL_KMIN)
            for kmax in sorted(data[a][V][l][kmin]):
              self._addRow(row, kmax, LVL_KMAX)
              for df in sorted(data[a][V][l][kmin][kmax]):
                self._addRow(row, df, LVL_DF)
                if not ignoreN:
                  Nstr = ""
                  for N in sorted(data[a][V][l][kmin][kmax][df]):
                    Nstr += str(N) + " "
                  self._addRow(row, Nstr, LVL_N)
                table.append(row)
                row = []
    print tabulate(table, headers="firstrow", floatfmt=".9f")
    
  def listData(self, data):
    for a in sorted(data):
      for V in sorted(data[a]):
        for l in sorted(data[a][V]):
          for kmin in sorted(data[a][V][l]):
            for kmax in sorted(data[a][V][l][kmin]):
              for df in sorted(data[a][V][l][kmin][kmax]):
                for N in sorted(data[a][V][l][kmin][kmax][df]):
                  print "\n" + makePathStr(a,V,l,kmin,kmax,df,N)
                  base = data[a][V][l][kmin][kmax][df][N]
                  if len(base[INDEX_ROOT]) > 0:
                    print "Roots:"
                    table = []
                    for value in base[INDEX_ROOT]:
                      table.append([self._getPrintVal(value.k), self._getPrintVal(value.E)])
                    print tabulate(table, floatfmt=".9f")
                  if len(base[INDEX_POLE]) > 0:
                    print "Poles:"
                    table = []
                    for value in base[INDEX_POLE]:
                      table.append([self._getPrintVal(value.k), self._getPrintVal(value.E)])
                    print tabulate(table, floatfmt=".9f")
  
  def listListData(self, listData):
    for searchVal in sorted(listData,cmp=lambda x,y: cmp(abs(x), abs(y))):
      print "\n**************************************"
      print "******* " + str(searchVal) + " *******"
      print "**************************************"
      self.listData(listData[searchVal])

###########
      
  def tabulateListSearch_k(self, data, lockLvl):
    self._tabulateListSearch(data, lockLvl, True)
    
  def tabulateListSearch_E(self, data, lockLvl):
    self._tabulateListSearch(data, lockLvl, False)

  def _tabulateListSearch(self, data, lockLvl, ks):
    table, colVals = self._createDataTable(data, lockLvl, ks)
    header = self._makeHeader(colVals, lockLvl)
    print tabulate([header]+table, headers="firstrow", floatfmt=".9f")

###########
    
  def plotListSearch_k(self, data, lockLvl, plotLvl):
    return self._plotListSearch(data, lockLvl, plotLvl, True)
    
  def plotListSearch_E(self, data, lockLvl, plotLvl):
    return self._plotListSearch(data, lockLvl, plotLvl, False)

  def _plotListSearch(self, data, lockLvl, plotLvl, ks):
    table, colVals = self._createDataTable(data, lockLvl, ks)
    xs = []
    ys = []
    for i in range(len(colVals)):
      xs.append(colVals[i][plotLvl])
    
    row = 0    
    for tabRow in table:
      ys.append([])
      for i in range(len(colVals)):
        col = i+1
        if tabRow[col] != "":
          ys[row].append(float(tabRow[col]))
        else:
          ys[row].append(np.nan)
      row += 1
    return (xs, ys) 
  
###########
  
  def totalListSearch_k(self, data, lockLvl, plotLvl):
    return self._totalListSearch(data, lockLvl, plotLvl, True)
    
  def totalListSearch_E(self, data, lockLvl, plotLvl):
    return self._totalListSearch(data, lockLvl, plotLvl, False)

  def _totalListSearch(self, data, lockLvl, plotLvl, ks):
    table, colVals = self._createDataTable(data, lockLvl, ks)
    xs = []
    ys = []
    for i in range(len(colVals)):
      xs.append(colVals[i][plotLvl])
      num = self._countNumInCol(table,i+1)
      ys.append(float(num))
    return (xs, ys)
      
###########
      
  def _countNumInCol(self, table, col):
    num = 0
    for row in table:
      if row[col] != "":
        num += 1
    return num

  def _createDataTable(self, data, lockLvl, ks):
    colMapper = {}
    colCnt = 1
    colVals = []
    table = []
    uniqueRowVals = self._getUniqueRowVals(data, lockLvl)
    for lockVal in sorted(uniqueRowVals):
      row = [""]*(colCnt)
      row[0] = str(lockVal)
      for searchVal in sorted(data,cmp=lambda x,y: cmp(abs(x), abs(y))):
        colCnt = self._R_tabulateListSearch(data[searchVal], row, lockLvl, lockVal, 0, [], ks, colMapper, colVals, colCnt)
      table.append(row)
    
    for row in table:
      while len(row)!=colCnt:
        row.append("")
    return (table, colVals)
  
  def _R_tabulateListSearch(self, data, row, lockLvl, lockVal, lvlCnt, lvlVals, ks, colMapper, colVals, colCnt):
    for lvlVal in sorted(data):
      if lockLvl!=lvlCnt or lockVal==lvlVal:
        nextlvlVals = lvlVals + [lvlVal]
        if lvlCnt != LVL_LAST:
          colCnt = self._R_tabulateListSearch(data[lvlVal], row, lockLvl, lockVal, lvlCnt+1, nextlvlVals, ks, colMapper, colVals, colCnt)
        else:
          pathStr = makePathStrFromLst(nextlvlVals,lockLvl)
          base = data[lvlVal]
          for typeIndex in range(len(base)):
            if len(base[typeIndex]) > 0:
              pathStr += "_" + getTypeStr(typeIndex)
              valStr = "" #Hopefully they'll only ever be one, if not just put all in the same cell.
              for value in base[typeIndex]:
                if ks:
                  valStr += self._getPrintVal(value.k) + " "
                else:
                  valStr += self._getPrintVal(value.E) + " "
              valStr = valStr.strip()    
              if pathStr in colMapper:
                row[colMapper[pathStr]] = valStr
              else:
                row.append(valStr)
                colVals.append(nextlvlVals)
                colMapper[pathStr] = colCnt #First index contains the searchVal
                colCnt += 1
    return colCnt
        
  def _getUniqueRowVals(self, data, lockLvl):
    uniqueVals = Set()
    for searchVal in data:
      for a in data[searchVal]:
        if lockLvl==LVL_A: uniqueVals.add(a)
        for V in data[searchVal][a]:
          if lockLvl==LVL_V: uniqueVals.add(V)
          for l in data[searchVal][a][V]:
            if lockLvl==LVL_L: uniqueVals.add(l)
            for kmin in data[searchVal][a][V][l]:
              if lockLvl==LVL_KMIN: uniqueVals.add(kmin)
              for kmax in data[searchVal][a][V][l][kmin]:
                if lockLvl==LVL_KMAX: uniqueVals.add(kmax)
                for df in data[searchVal][a][V][l][kmin][kmax]:
                  if lockLvl==LVL_DF: uniqueVals.add(df)
                  for N in data[searchVal][a][V][l][kmin][kmax][df]:
                    if lockLvl==LVL_N: uniqueVals.add(N)
    return list(uniqueVals)
  
  def _makeHeader(self, colVals, lockLvl):
    header = ["x"]
    for colVal in colVals:
      header.append(makePathStrFromLst(colVal,lockLvl))
    return header

  def _addRow(self, row, val, num):
    if len(row)==0:
      for i in range(num):
        row.append("")
    row.append(str(val))

  def _getPrintVal(self, val):
    if self.printReal and self.printImag:
      ret = str(val)
    elif self.printReal:
      ret = '%.9f' % val.real
    else:
      ret = '%.9f' % val.imag
    return ret
  
###########
  
  def tabulatePoleData(self,data,a,V,l,kmin,kmax,df,forExcel=True,useEnergy=False):
    t = Tables(useEnergy)
    base = data[a][V][l][kmin][kmax][df]
    print "a="+str(a)+"_V="+str(V)+"_l="+str(l)+"_kmin="+str(kmin)+"_kmax="+str(kmax)+"_df="+str(df)
    t.tabulatePoleData(base,forExcel)
              
###########

  def plotPoles_E(self,data,a=None,V=None,l=None,kmin=None,kmax=None,df=None,N=None):
    xs = []
    ys = []
    if a is None:
      base = data
      for val in base:
        for pole in base[val][V][l][kmin][kmax][df][N][INDEX_POLE]:
          self._addVals(xs,ys,val,pole)
    elif V is None:
      base = data[a]
      for val in base:
        for pole in base[val][l][kmin][kmax][df][N][INDEX_POLE]:
          self._addVals(xs,ys,val,pole)
    elif l is None:
      base = data[a][V]
      for val in base:
        for pole in base[val][kmin][kmax][df][N][INDEX_POLE]:
          self._addVals(xs,ys,val,pole)
    elif kmin is None:
      base = data[a][V][l]
      for val in base:
        for pole in base[val][kmax][df][N][INDEX_POLE]:
          self._addVals(xs,ys,val,pole)
    elif kmax is None:
      base = data[a][V][l][kmin]
      for val in base:
        for pole in base[val][df][N][INDEX_POLE]:
          self._addVals(xs,ys,val,pole)
    elif df is None:
      base = data[a][V][l][kmin][kmax]
      for val in base:
        for pole in base[val][N][INDEX_POLE]:
          self._addVals(xs,ys,val,pole)
    elif N is None:
      base = data[a][V][l][kmin][kmax][df]
      for val in base:
        for pole in base[val][INDEX_POLE]:
          self._addVals(xs,ys,val,pole)
    return (xs, ys)
          
  def _addVals(self,xs,ys,val,pole):
    xs.append(float(val))
    ys.append(pole.E.real)
      
def getData(dataDir):
  d = Dat()
  print "\nGetting Data from dir" + dataDir
  for root, dirs, files in os.walk(dataDir, topdown=False):
    for fileName in files:
        try:
            name = os.path.splitext(fileName)[0]
            if "Coeff" not in name:
              path = [float(x) if "." in x else int(x) for x in name.split("_")[1:]]
              p = Parser()
              with open(dataDir+fileName, 'r') as f:
                for line in f:
                  #print line
                  if p.isN(line):
                    N = p.getN(line)
                  elif p.isRoot(line):
                    p.setRoot(d, line, path, N)
                  elif p.isPole(line):
                    p.setPole(d, line, path, N)  
        except Exception as e:
            print "\tError Parsing: " + fileName + "  " + str(e)
  print "\tFinished"
  return d
  

FILE_STR = base+"./Data/0.01/1.0/"
if __name__ == "__main__":
  d = getData(FILE_STR)
  p = Printer(True, False)
  f = Filter()
  print d.getRange(LVL_KMIN)
  p.listParameters(d.getData())
  p.listData(d.getData())
  p.listData(d.searchkForRoot(-5.48365-1.41865j))
  p.listData(d.searchEForRoot(14.02891+7.77940j))
  p.listData(d.searchkForPole(-0.0-0.60917j))
  p.listData(d.searchEForPole(1.57091+2.12479j))
  p.listListData(d.listSearch(d.searchEForPole,[-0.102053+0j, -0.162536+0j]))
  f.selectWhereIn(d.getData(),N_s=[65])
  p.tabulateListSearch_E(d.listSearch(d.searchEForPole,[-0.102053+0j, -0.162536+0j]),LVL_V)