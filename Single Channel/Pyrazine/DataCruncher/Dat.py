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
sys.path.insert(0,base+'/..')
import General.Numerical as num
from GenDat import *
from PyrazineDataVals import *

INDEX_EMIN = 0
INDEX_EMAX = 1
INDEX_OFFSET = 2
INDEX_DF = 3
INDEX_ERANGE = 4

################### Init Functions ###################
    
class Dat:
  def __init__(self):
    self.data = {}

################### Data Build Functions ###################
    
  def addRoot(self,key,N,k,E,prevRootVal=None,diff=None):
    end = self._addPathToSelfData(key,N)
    end[INDEX_ROOT].append(Root(k,E,prevRootVal,diff))
  
  def addPole(self,key,N,k,E,isNew,isLost):
    end = self._addPathToSelfData(key,N)
    end[INDEX_POLE].append(Pole(k,E,isNew,isLost))
  
  def _addPathToSelfData(self,key,N):
    return self._addPath(self.data,key,N)
  
  def _addPath(self,data,key,N):
    if key not in data:
      data[key] = {}
    if N not in data[key]:
      data[key][N] = [[],[],[]]
    end = data[key][N]
    return end
  
####################### Accessor #######################
  
  def getData(self):
    return self.data

###########

class Printer:  
  def tabulatePoleData(self,data,dataVal,forExcel=False,useEnergy=True):
    t = Tables(useEnergy)
    base = data[dataVal]
    t.tabulatePoleData(base,forExcel)
              
###########
  
  
class Parser:                  
  def isN(self, line):
    return line[0] == 'N'

  def getN(self, line):
    return int(line[2:line.find(',')])

  def isRoot(self, line):
    return "Root_" in line
    
  def setRoot(self, d, line, key, N):
    k,E = self._extractNums(line)
    prevRootVal, diff = self._extractRootInfo(line)
    d.addRoot(key,N,k,E,prevRootVal,diff)

  def isPole(self, line):
    return "Pole_" in line and "MiscPole_" not in line

  def setPole(self, d, line, key,N):
    new = POLE_STATUS_NEW in line
    lost = POLE_STATUS_LOST in line
    k,E = self._extractNums(line)
    d.addPole(key,N,k,E,new,lost)
    
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

def _getKey(path, pathIndex):
    if pathIndex == INDEX_ERANGE:
        return path[INDEX_EMAX] - path[INDEX_EMIN]
    else:
        return path[pathIndex]
  
def getData(dataDir, pathIndex):
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
                    p.setRoot(d, line, _getKey(path,pathIndex), N)
                  elif p.isPole(line):
                    p.setPole(d, line, _getKey(path,pathIndex), N)  
        except Exception as e:
            print "\tError Parsing: " + fileName + "  " + str(e)
  print "\tFinished"
  return d
