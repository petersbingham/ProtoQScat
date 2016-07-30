from ratsmatwrap import *
from general.numerical import *
import argparse

args_ephasefit = argparse.ArgumentParser(description="Eigenphase fit", add_help=False)
args_ephasefit.add_argument("plotStart_", help="Plot Start Energy", type=float, default=None)
args_ephasefit.add_argument("plotEnd_", help="Plot End Energy", type=float, default=None)
args_ephasefit.add_argument("steps_", help="Number of steps", type=int, default=None)
args_ephasefit.add_argument("polyOrder_", help="Order of Poly Fit", type=int, default=2)
args_ephasefit.add_argument("resDatas_", help="Poles to fit to", type=complex, nargs="+")
args = args_ephasefit.parse_args()

FIT_CYCLE = ['red', 'blue', 'purple', 'orange', 'cyan']

try:
    kmat = RatSMatWrap(FILENAME)
    ratEPhaseMats, ePhaseFitPointsLst = kmat.getEigenSumFit(args.resDatas_, args.polyOrder_, args.plotStart_, args.plotEnd_, args.steps_)
    
    ratEPhaseMats.setDetails("R-matrix", ['green'])
    for i in range(len(ePhaseFitPointsLst)):
        ePhaseFitPointsLst[i].setDetails(str(truncateComplex(6,args.resDatas_[i])), FIT_CYCLE[i])
    
    sm.plotSum("Lin", [ratEPhaseMats]+ePhaseFitPointsLst, title=DESC_STR+" Fitted Resonance Model", xlabel="Electron Energy (rydbergs)", ylabel="Total Phase") 
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()