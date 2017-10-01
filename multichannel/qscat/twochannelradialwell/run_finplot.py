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
spArgs.add_argument("m_", help="Element m", type=int)
spArgs.add_argument("n_", help="Element n", type=int)
args = spArgs.parse_args()

try:
    kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, massMult=MASSMULT)
    matSeq = sm.matSequence() 
    anaSmat, ratSmat = getSmats(args, kCal, kCal)
    
    vals = ratSmat.getFinRange(args.startE_, args.endE_, args.plotComplex_, args.steps_, args.m_, args.n_)
    
    sp.plotSingle("Two Channel Radial Well Fit : Fin", vals[0], [vals[1],vals[2]], 'Total Energy (hartrees)', 'Fin', ['real','imag'])

except (DCException, sm.MatException) as inst:
    print str(inst)
    sys.exit()