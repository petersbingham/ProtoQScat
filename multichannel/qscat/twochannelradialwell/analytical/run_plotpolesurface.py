from runargscommon import *
from general import *
import general.numerical as num
from tabulate import tabulate
import numpy as np

rootArgs = argparse.ArgumentParser(description="Two Channel Radial Well - Find Resonant Roots", parents=[parentArgs])
args = rootArgs.parse_args()

kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, eneFactor=ENEFACTOR)  
mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
smat = Smat(args.r0_, mats)

Rs = np.arange(-10, 200, 0.2) #0.1, 0.2, 0.4
Is = np.arange(-100, 100, 0.4)
smat.plotPoles(Rs, Is, True)