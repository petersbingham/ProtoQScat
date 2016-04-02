from Analytical.runArgsRange import *
from runBase import *
from cmath import *

spArgs = argparse.ArgumentParser(description="Two Channel Radial Well Plot", parents=[tcp_range])
spArgs.add_argument("plotStart_", help="Plot Start", type=float)
spArgs.add_argument("plotEnd_", help="Plot End", type=float)
spArgs.add_argument("plotComplex_", help="Plot Complex Energy Offset", type=float)
spArgs.add_argument("plotSteps_", help="Plot Steps", type=int)
spArgs.add_argument("plotType_", help="Plot Type")
args = spArgs.parse_args()

def initPlotVars():
    if args.plotType_=="Lin":
      x = args.plotStart_+args.plotComplex_*1.0j
      dx = (args.plotEnd_-args.plotStart_) / args.plotSteps_
    else:
      x = log10(args.plotStart_)
      dx = (log10(args.plotEnd_)-log10(args.plotStart_)) / args.plotSteps_
    return x, dx
x, dx = initPlotVars()

def getParameterDesc(args):
    c = num.Compare()
    if not c.floatCompare(args.lam_,0.0):
      if not c.floatCompare(args.t1_,0.0) or not c.floatCompare(args.t2_,0.0):
        s = "Coupled Inelastic"
      else:
        s = "Coupled Elastic"
    else:
      if not c.floatCompare(args.t1_,0.0) or not c.floatCompare(args.t2_,0.0):
        s = "Uncoupled Inelastic"
      else:
        s = "Uncoupled Elastic"
    return s + " (" + str([args.r0_, args.v1_, args.v2_, args.t1_, args.t2_, args.lam_]) + ") Radial Well "

def doMatPlot(args, mat, imag, title, signString):
    matSeq = sm.matSequence()

    x, dx = initPlotVars()
    for i in range(0,args.plotSteps_+1,1):
      ene = getEnergy(args.plotType_, x)
      mat.setEnergy(ene)
      matSeq[ene] = mat.getMatrix()
      x += dx
      
    matSeq.setDetails("", ['red', 'green', 'blue', 'purple'])

    sm.plotAll(args.plotType_, [matSeq], imag, title+"\n"+signString)