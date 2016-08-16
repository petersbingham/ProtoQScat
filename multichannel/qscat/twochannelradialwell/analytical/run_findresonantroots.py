from runargscommon import *
from general import *
import general.numerical as num
from tabulate import tabulate

rootArgs = argparse.ArgumentParser(description="Two Channel Radial Well - Find Resonant Roots", parents=[parentArgs])
args = rootArgs.parse_args()

kCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_COMP, eneFactor=ENEFACTOR)  
mats = Mats(args.v1_, args.v2_, args.lam_, kCal)
smat = Smat(args.r0_, mats)
    
rootStarts = [6.1314979779987873+5.9402238355294239j,
6.1314979779987873-5.9402238355294239j,
6.4820184783656583+7.1297150595251084j,
6.4820184783656583-7.1297150595251084j,
24.4816062940894668+14.6503227997078671j,
24.4816062940894668-14.6503227997078671j,
24.6968714409361247+16.6233906322054104j,
24.6968714409361247-16.6233906322054104j,
53.0405510180163589+24.5715845570572036j,
53.0405510180163589-24.5715845570572036j,
53.1688166730558152+27.3473846710057344j,
53.1688166730558152-27.3473846710057344j,
91.6626138248640956+35.3717803242041029j,
91.6626138248640956-35.3717803242041029j,
91.726360915098283+38.9523092399742978j,
91.726360915098283-38.9523092399742978j]

c = num.RationalCompare(1.0E-07, 0.1)

table = []
offsetInc = 0.001  
for rootStart in rootStarts:
    rootStartMod = num.truncateComplex(3, rootStart)
    found = False
    for offSetIndex in range(100):
        if found:
            break
        for dir in range(2):
            
            root = "none" 
            if dir==0:
                rootStartMod = rootStartMod.real-offSetIndex*offsetInc + rootStartMod.imag*1.0j
            else:
                rootStartMod = rootStartMod.real+offSetIndex*offsetInc + rootStartMod.imag*1.0j
                
            try:
                root = smat.findRoot(rootStartMod)
                if not c.isClose(root, rootStart):
                    root = "none"
                    continue
            except ValueError:
                continue
            found = True
            break
    print "-"    
    table.append([rootStartMod,rootStart,root])
    
print getFormattedHTMLTable(table, ["Search Starting Value","Numerical Pole","Search Return"], 9, numalign="decimal")