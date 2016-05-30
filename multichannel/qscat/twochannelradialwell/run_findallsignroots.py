from runbase import *
from analytical.runargsrange import *
import scattering.stran as S
import general.numerical as num
import tabulate

spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Fit - Pole find", parents=[tcp_range])
spArgs.add_argument("startE_", help="Start Energy", type=complex)
args = spArgs.parse_args()

anaSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]
fitSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]
ratSignList = [[1.0,1.0],[1.0,-1.0],[-1.0,1.0],[-1.0,-1.0]]

table = [["Ana Signs","Fit Signs","Jost Signs","1st Starting Point","2nd Starting Point"]]

def getRootStr(root):
    if root is not None:
        return "{:.8f}".format(root)
    else:
        return "none\t"

def getRoot(ratkCal, ratSmat, signString):
    row = signString
    #row = [" ", "|", " "]
    roots = ratSmat.findConjRoots(args.startE_)
    row.append(getRootStr(roots[0]))
    cmp = num.Compare()
    if roots[1] is not None and roots[0] is not None and cmp.complexCompare(roots[1], roots[0]):
        row.append("same")
    else:
        row.append(getRootStr(roots[1]))
    table.append(row)
   
dokSignIt(args, anaSignList, fitSignList, ratSignList, None, getRoot, False, True)
print tabulate.tabulate(table,headers="firstrow", tablefmt="html")