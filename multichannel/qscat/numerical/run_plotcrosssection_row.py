from ratsmatwrap import *
import matplotlib.pyplot as plt
from runbase import *

childArgs = argparse.ArgumentParser(parents=[parentArgs])
childArgs.add_argument("m_", help="m element", type=int)
args = childArgs.parse_args()

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  xs = kmat.getDiscreteXS()
  ratxs = kmat.getRatXS()
  
  xs.setDetails("Raw", ['orange','brown','cyan'])
  ratxs.setDetails("Ana", ['red','green','blue'])
  
  sm.plotRow("Lin", [xs, ratxs], args.m_, False, "Pyrazine Cross Sections", 'Electron Energy (rydbergs)', 'Cross Section (bohr^2)')
  
  plt.show()
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()