import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')
from RatSMat import *
from Analytical.DoubleChannel import *
from PoleFinder import *

def _getTypeName(args, anakCal, fitkCal):
    return "Two Channel Radial Well_" + str(anakCal) + "_" + str(fitkCal) + "_" + str(args.r0_) + "_" + str(args.v1_) + "_" + str(args.v2_) + "_" + str(args.lam_) + "_" + str(args.eneStart_) + "_" + str(args.eneEnd_) + "_" + str(args.eneComplex_) + "_" + str(args.eneSteps_)

def getEnergy(type, x): 
    if type == "Lin":
        return x
    elif type=="Log" or type=="LogLog":
        return pow(10.0, x)

def getAnaSmat(args, anakCal):  
    mats = Mats(args.v1_, args.v2_, args.lam_, anakCal)
    return Smat(args.r0_, mats)

def getDiscreteAnaSmats(args, anaSmat):
    dEne = (args.eneEnd_-args.eneStart_) / float(args.eneSteps_)
    ene = args.eneStart_ + args.eneComplex_*1.0j
    sMats = {}
    for i in range(0,args.eneSteps_+1,1):
        anaSmat.setEnergy(ene)
        sMats[ene] = anaSmat.getMatrix()
        ene += dEne
    return sMats 

def getRatSmat(args, anaSmat, anakCal, fitkCal, suppressCmdOut=False):
    sMats = getDiscreteAnaSmats(args, anaSmat)  
    return RatSMat(sMats, fitkCal, fitName=_getTypeName(args, anakCal, fitkCal), suppressCmdOut=suppressCmdOut)

def getSmats(args, anakCal, fitkCal):
    anaSmat = getAnaSmat(args, anakCal)
    ratSmat = getRatSmat(args, anaSmat, anakCal, fitkCal)
    return (anaSmat, ratSmat)

def getDecimatedRatSmat(args, smats, anakCal, fitkCal, N, suppressCmdOut=False): #To get ratSmat for same data set as for getPolyRoots to allow comparison with Muller
    decimator = Decimator(0, len(smats)-1, 0)
    sMats, decStr = decimator.decimate(smats, N)
    return RatSMat(sMats, fitkCal, fitName=_getTypeName(args, anakCal, fitkCal), suppressCmdOut=suppressCmdOut)

def getPolyRoots(args, anakCal, fitkCal, resultsPath, mode, cmpValue=None):
    anaSmat = getAnaSmat(args, anakCal)
    smats = getDiscreteAnaSmats(args, anaSmat)
    PoleFinder(smats, fitkCal, resultsPath, _getTypeName(args, anakCal, fitkCal), ENEFACTOR, 0, len(smats)-1, 0, 0.05, 1, cmpValue=cmpValue, mode=mode)
      
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
                        ratSmat.kCal = ratkCal
                        if signsAsList:
                            signs = [str(anakCal), str(fitkCal), str(ratkCal)]
                        else:
                            signs = "ana:"+str(anakCal)+", fit:"+str(fitkCal)+", jost:"+str(ratkCal)
                        ratFun(ratkCal, ratSmat, signs)
    except (DCException, sm.MatException) as inst:
        print str(inst)
        sys.exit()
