from RatSMatWrap import *
import matplotlib.pyplot as plt
from runbase import *
args = parentArgs.parse_args()

sm.setParams(0.12, 0.07, 0.90, 0.93, 0.20, 0.11)
sm.setSubPlots(2,1,"Pyrazine Cross Sections", 'Electron Energy (rydbergs)', 'Cross Section (bohr^2)')

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  totXss = kmat.getTotalDiscreteXS()
  totXss.setConversionEnergy(sm.eVs)
  totXss.writeAllToFile()
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()