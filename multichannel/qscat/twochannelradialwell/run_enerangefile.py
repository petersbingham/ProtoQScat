from runbase import *
from analytical.runargsrange import *
import scattering.stran as S

spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - S-matrix Range", parents=[tcp_range])
spArgs.add_argument("rangeStart_", help="Range Start", type=float)
spArgs.add_argument("rangeEnd_", help="Range End", type=float)
spArgs.add_argument("rangeComplex_", help="Range Complex Offset", type=float)
spArgs.add_argument("rangeSteps_", help="Range Steps", type=int)
spArgs.add_argument("rangeType_", help="Range Type")
args = spArgs.parse_args()

x = args.rangeStart_+args.rangeComplex_*1j
dx = (args.rangeEnd_-args.rangeStart_) / args.rangeSteps_

try:
    kCal = sm.kCalculator([args.t1_,args.t2_], eneFactor=ENEFACTOR)
    matSeq = sm.matSequence() 
    anaSmat, ratSmat = getSmats(args, kCal, kCal)
    
    if args.rangeType_ == "S":
        rangeMat = ratSmat
    elif args.rangeType_ == "UnitaryOp": 
        rangeMat = S.UniOpmat(ratSmat)
      
    for i in range(0,args.rangeSteps_+1,1):
        ene = getEnergy("Lin", x)
        rangeMat.setEnergy(ene)
        matSeq[ene] = rangeMat.getMatrix()
        x += dx
    matSeq.writeAllToFile()
  
except (DCException, sm.MatException) as inst:
    print str(inst)
    sys.exit()