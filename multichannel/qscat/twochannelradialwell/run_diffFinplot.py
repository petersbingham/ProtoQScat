from runbase import *
import matplotlib.pyplot as plt
from analytical.runargsrange import *
from scattering.stran import *
import general.simpplot as sp

spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Plot det(Fin)", parents=[tcp_range])
spArgs.add_argument("startE_", help="Start Energy", type=float)
spArgs.add_argument("endE_", help="End Energy", type=float)
spArgs.add_argument("plotComplex_", help="Plot Complex Energy Offset", type=float)
spArgs.add_argument("steps_", help="Number of steps", type=int)
spArgs.add_argument("m_", help="Matrix m", type=int)
spArgs.add_argument("n_", help="Matrix n", type=int)
args = spArgs.parse_args()

try:
    kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, massMult=MASSMULT)
    matSeq = sm.matSequence() 
    anaSmat, ratSmat = getSmats(args, kCal, kCal)
    
    vals = ratSmat.getFinRange(args.startE_, args.endE_, args.plotComplex_, args.steps_, args.m_, args.n_)
    dvals = ratSmat.getDiffFinRange(args.startE_, args.endE_, args.plotComplex_, args.steps_, args.m_, args.n_)
    
    sp.plotSingle("Two Channel Radial Well Fit : real Fin", vals[0], [vals[1],dvals[1]], 'Total Energy (hartrees)', 'Fin', ['function','differential'],markerSz=None)
    sp.plotSingle("Two Channel Radial Well Fit : imag Fin", vals[0], [vals[2],dvals[2]], 'Total Energy (hartrees)', 'Fin', ['function','differential'],markerSz=None)

except (DCException, sm.MatException) as inst:
    print str(inst)
    sys.exit()