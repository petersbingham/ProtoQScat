from ratsmatwrap import *
import matplotlib.pyplot as plt
import argparse
parentArgs = argparse.ArgumentParser(description="Pyrazine Fit Routine", add_help=False)
parentArgs.add_argument("subStart_", help="Sub Set Start", type=int, nargs='?', default=-1)
parentArgs.add_argument("subEnd_", help="Sub Set End", type=int, nargs='?', default=-1)
args = parentArgs.parse_args()

KNOWN = True

NUM_FND = 0

def _findRoot(N, starting):
  kmat = RatSMatWrap("fort.19",N,args.subStart_,args.subEnd_, kfitSigns=[1.0,1.0,1.0])
  new = round(starting.real,3) + round(starting.imag,4)*1.0j
  root = kmat.findRoot_(new)
  print str(new) + "   " + str(root) + "\n\n"
          
def _findKnownRoot(N):
  kmat = RatSMatWrap("fort.19",N,args.subStart_,args.subEnd_, kfitSigns=[1.0,1.0,1.0])
  starting = 0.076641273-0.0006538903j
  print str(starting) + "   " + str(kmat.findRoot(starting)) + "\n\n"
          
try:
    if KNOWN:
        _findKnownRoot(4)
        _findKnownRoot(6)
        _findKnownRoot(8)
        _findKnownRoot(10)
        _findKnownRoot(12)
        _findKnownRoot(14)
        _findKnownRoot(16)
        _findKnownRoot(18)
        _findKnownRoot(20)
        _findKnownRoot(22)
        _findKnownRoot(24)
        _findKnownRoot(26)
        _findKnownRoot(28)
        _findKnownRoot(30)
    else:
        _findRoot(4, 0.08508759093730-0.00361311967739j)
        _findRoot(6, 0.07663626762359-0.00064893798744j)
        _findRoot(8, 0.07736096864253-0.00064262088651j)
        _findRoot(10, 0.07536112329385-0.00040573838314j)
        _findRoot(12, 0.07865655198363-0.00181817072217j)
        _findRoot(14, 0.04736396051274-0.01199006660536j)
        _findRoot(16, 0.08920346198933-0.03314890230006j)
        _findRoot(18, 0.08246371483663-0.01179852552228j)
        _findRoot(20, 0.06184249601552+0.00236666865138j)
        _findRoot(22, 0.06601948460344-0.02098239655253j)
        _findRoot(24, 0.01425703435212+0.00820492609632j)
        _findRoot(26, 0.01515240507248+0.00103440892864j)
        _findRoot(28, 0.00162780117706+0.00153073193941j)
        _findRoot(30, 0.02350318740315-0.00755343673500j)
except (sm.MatException) as inst:
  print str(inst)
  sys.exit()