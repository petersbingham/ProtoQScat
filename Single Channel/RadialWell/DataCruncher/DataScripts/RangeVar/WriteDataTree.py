import sys
import os
rangeDir = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,rangeDir+'./..')
from Base import *

dataDir = rangeDir+'/'+sys.argv[1]+'/'
print dataDir
d = getData(dataDir)
p = Printer(True, False)
p.listParameters(d.getData())