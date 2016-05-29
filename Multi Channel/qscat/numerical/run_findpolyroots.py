from RatSMatWrap import *
import matplotlib.pyplot as plt
from runbase import *
args = parentArgs.parse_args()

try:
  kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
  roots = kmat.findPolyRoots()
  for root in roots:
      print root
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()