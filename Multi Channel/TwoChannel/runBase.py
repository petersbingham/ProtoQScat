import sys
#sys.path.append("../../Utilities")
sys.path.append("..")
from RatSMat import *
from Analytical.DoubleChannel import *

def _getTypeName(anakCal, fitkCal, r0, v1, v2, lam, eneStart, eneEnd, eneComplex, eneSteps):
  return "_" + str(anakCal) + "_" + str(fitkCal) + "_" + str(r0) + "_" + str(v1) + "_" + str(v2) + "_" + str(lam) + "_" + str(eneStart) + "_" + str(eneEnd) + "_" + str(eneComplex) + "_" + str(eneSteps)

def getEnergy(type, x): 
  if type == "Lin":
    return x
  elif type=="Log" or type=="LogLog":
    return pow(10.0, x)

def getAnaSmat(args, anakCal):  
  mats = Mats(args.v1_, args.v2_, args.lam_, anakCal)
  return Smat(args.r0_, mats)

def getRatSmat(args, anaSmat, anakCal, fitkCal, suppressCmdOut=False):
  dEne = (args.eneEnd_-args.eneStart_) / float(args.eneSteps_)
  ene = args.eneStart_ + args.eneComplex_*1.0j
  sMats = {}
  for i in range(0,args.eneSteps_+1,1):
    anaSmat.setEnergy(ene)
    sMats[ene] = anaSmat.getMatrix()
    ene += dEne  
  return RatSMat(sMats, fitkCal.kl, fitName="Two Channel Radial Well" + _getTypeName(anakCal, fitkCal, args.r0_, args.v1_, args.v2_, args.lam_, args.eneStart_, args.eneEnd_, args.eneComplex_, args.eneSteps_), suppressCmdOut=suppressCmdOut)

def getSmats(args, anakCal, fitkCal):
  anaSmat = getAnaSmat(args, anakCal)
  ratSmat = getRatSmat(args, anaSmat, anakCal, fitkCal)
  return (anaSmat, ratSmat)
  
def dokSignIt(args, anaSignList, fitSignList, ratSignList, anaFun, ratFun, suppressCmdOut=False, signsAsList=False):
    try:
        for anaSigns in anaSignList:
            anakCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=anaSigns, eneFactor=ENEFACTOR)
            anaSmat = getAnaSmat(args, anakCal)
            if signsAsList:
                signs = [str(anakCal)]
            else:
                signs = "ana:"+str(anakCal)
            if anaFun is not None:
                anaFun(anakCal, anaSmat, signs)
            if ratFun is not None:
                for fitSigns in fitSignList:
                    fitkCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=fitSigns, eneFactor=ENEFACTOR)
                    ratSmat = getRatSmat(args, anaSmat, anakCal, fitkCal, suppressCmdOut)
                    for ratSigns in ratSignList:
                        ratkCal = sm.kCalculator([args.t1_,args.t2_], ktype=sm.K_SIGN, ksigns=ratSigns, eneFactor=ENEFACTOR)
                        ratSmat.kFun = ratkCal.kl
                        if signsAsList:
                            signs = [str(anakCal), str(fitkCal), str(ratkCal)]
                        else:
                            signs = "ana:"+str(anakCal)+", fit:"+str(fitkCal)+", jost:"+str(ratkCal)
                        ratFun(ratkCal, ratSmat, signs)
    except (DCException, sm.MatException) as inst:
        print str(inst)
        sys.exit()
