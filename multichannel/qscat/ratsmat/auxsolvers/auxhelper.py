import numpy as np
import scipy.sparse.linalg as sp_sparse_linalg
import sympy as sy
from sympy.matrices import Matrix as sy_matrix
import sympy.polys as sy_polys
from sympy import poly
import mpmath
import collections

import sys
import os
base =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../../utilities')
sys.path.insert(0,base+'/..')

from general import *
import general.type_wrap as tw
import general.numerical as num
import scattering.matrices as sm

def canCacheCoefficients():
    if tw.mode == tw.mode_norm:
        npVer = [int(val) for val in np.__version__.split('.')]
        if npVer[0]>1 or (npVer[0]==1 and npVer[1]>6) or (npVer[0]==1 and npVer[1]==6 and npVer[2]>1): #saveTxt not supported prior.
            return True
    else:
        return True
    return False

class AuxHelper():
    def __init__(self, suppressCmdOut):
        self.printed = False
        self.suppressCmdOut = suppressCmdOut
    
    def setResultFileHandler(self, resultFileHandler):
        self.resultFileHandler = resultFileHandler
    
    def _printCalStr(self, nounStr, verbStr, wereLoaded):
        if not self.suppressCmdOut and not self.printed:
            addStr = ""
            if wereLoaded:
                addStr = "had been "
            print nounStr + " " + addStr + verbStr + " using " + self.typeStr 
            self.printed = True

    def _startLogAction(self, string):
        if self.resultFileHandler:
            self.resultFileHandler.startLogAction(string)

    def _endLogAction(self, string):
        if self.resultFileHandler:
            self.resultFileHandler.endLogAction(string)