from tabulate import tabulate
import General.Numerical as num
INDEX_ROOT = 0
INDEX_POLE = 1

Ns = [5,9,17,33,65]

def getTypeStr(index):
  if index==INDEX_ROOT:
    return "Root"
  if index==INDEX_POLE:
    return "Pole"

POLE_STATUS_NEW = "NEW"
POLE_STATUS_LOST = "LOST"
POLE_STATUS_REP = "REP"

class Val:
  def __init__(self, k, E):
    self.k = k
    self.E = E
class Root(Val):
  def __init__(self, k, E, prevRootVal, diff):
    Val.__init__(self, k, E)
    if prevRootVal:
      self.prevRootInd = prevRootVal-1
    else:
      self.prevRootInd = None
    if diff is not None:
      self.diff = diff
    else:
      self.diff = None
class Pole(Val):
  def __init__(self, k, E, isNew, isLost):
    Val.__init__(self, k, E)
    self.isNew = isNew
    self.isLost = isLost
  def getStatus(self):
    if self.isNew:
      return POLE_STATUS_NEW
    elif self.isLost:
      return POLE_STATUS_LOST
    else:
      return POLE_STATUS_REP


class Tables:
    def __init__(self, useEnergy=False):
      self.useEnergy = useEnergy
        
    def tabulatePoleData(self,base,forExcel=True):
      Ntable = [["N","status","(N)Root.r","(N)Root.i","(N-x)Root.r","(N-x)Root.i","diff.r","diff.i"]]
      Ptable = [["N","status","(N)Root.r","(N)Root.i","(N-x)Root.r","(N-x)Root.i","diff.r","diff.i"]]
      polesList = []
      prevNRoots = None
      for N in sorted(base):
        roots = base[N][INDEX_ROOT]
        if prevNRoots is not None:
          if len(base[N][INDEX_POLE]) > 0:
            first = True
            poleCnt = 0
            for pole in base[N][INDEX_POLE]:
              row = []
              if first:
                row.append(str(N))
                first = False
              else:
                row.append("")
              status = pole.getStatus()
              row.append(status)
              if status==POLE_STATUS_NEW or status==POLE_STATUS_REP:
                self._appendComplex(row, pole)
                root, prevRoot = self._getRoots(pole, roots, prevNRoots)
                if root is not None:
                    self._appendComplex(row, prevRoot)
                    self._appendComplex2(row, root.diff)
                else:
                    row.append("Conflict")
                    row.append("Conflict")
                    row.append("Conflict")
                    row.append("Conflict")
              elif status==POLE_STATUS_LOST:
                closestRoot = self._findClosestRoot(pole, roots)
                if closestRoot is not None:
                    self._appendComplex(row, closestRoot)
                else:
                    row.append("Conflict")
                    row.append("Conflict")
                self._appendComplex(row, pole)
                if not self.useEnergy:
                    self._appendComplex2(row, num.complexDiff(closestRoot.k,pole.k))
                else:
                    self._appendComplex2(row, num.complexDiff(closestRoot.E,pole.E))
              Ntable.append(row)
              row[0] = str(N)
              if poleCnt==len(polesList):
                polesList.append([])
              polesList[poleCnt].append(row)
              poleCnt += 1
          else:
            row = [str(N),"","","","","","",""]
            Ntable.append(row)
          self._insertTableDivide(Ntable, forExcel)
        prevNRoots = roots
      self._insertUnderLine(forExcel)
      print "\nGrouped_by_N:"
      self._tabulate(Ntable, forExcel)
      first = True
      for poleList in polesList:
        if not first:
          self._insertTableDivide(Ptable, forExcel)
        first = False
        for poleRow in poleList:
          Ptable.append(poleRow)
      print "\nGrouped_by_Pole:"
      self._tabulate(Ptable, forExcel)
    
    def _appendComplex(self, row, c):
      if not self.useEnergy:
        row.append(c.k.real)
        row.append(c.k.imag)
      else:
        row.append(c.E.real)
        row.append(c.E.imag)
    
    def _appendComplex2(self, row, c):
      row.append(c.real)
      row.append(c.imag)
        
    def _getRoots(self, pole, roots, prevNRoots):
      foundRoot = None
      for root in roots:
        if not self.useEnergy:
          if root.k == pole.k:
            foundRoot = root
        else:
          if root.E == pole.E:
            foundRoot = root
      if foundRoot is not None:
          return foundRoot, prevNRoots[foundRoot.prevRootInd]
      else:
          return None, None
    
    def _findClosestRoot(self, pole, roots):
      closestRoot = None
      for root in roots:
        if closestRoot is None:
          closestRoot = root  
        else:
          if not self.useEnergy:
            if num.absDiff(pole.k,root.k) < num.absDiff(pole.k,closestRoot.k):
              closestRoot = root
          else:
            if num.absDiff(pole.k,root.k) < num.absDiff(pole.k,closestRoot.k):
              closestRoot = root
      return closestRoot
    
    def _insertUnderLine(self, forExcel):
      if not forExcel:
          print "_______________________________________________________________________________________________"
    
    def _insertTableDivide(self, table, forExcel):
      table.append(["---","--------","---------------","---------------","---------------","---------------","------------","------------"])      
    
    def _tabulate(self, table, forExcel):
      if not forExcel: 
          table.append(["","","","","","","",""])
          print tabulate(table, headers="firstrow", floatfmt=".9f", tablefmt="html")
      else:
          print tabulate(table, headers="firstrow", floatfmt=".9f")