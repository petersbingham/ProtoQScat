import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/..')
from RatSMat import *
from Analytical.DoubleChannel import *
from PoleFinder import *
from ResultFileHandler import *

def _getFileHandler(args, anakCal, fitkCal):
    sysName = "Two Channel Radial Well_" + str(anakCal) + "_" + str(fitkCal) + "_" + str(args.r0_) + "_" + str(args.v1_) + "_" + str(args.v2_) + "_" + str(args.t1_) + "_" + str(args.t2_) + "_" + str(args.lam_) + "_" + str(args.eneStart_) + "_" + str(args.eneEnd_) + "_" + str(args.eneComplex_) + "_" + str(args.eneSteps_)
    return ResultFileHandler(sysName)

def getEnergy(type, x): 
    if type == "Lin":
        return x
    elif type=="Log" or type=="LogLog":
        return pow(10.0, x)

def getAnaSmat(args, anakCal):  
    mats = Mats(args.v1_, args.v2_, args.lam_, anakCal)
    return Smat(args.r0_, mats)

def getDiscreteAnaSmats(args, anaSmat=None):
    dEne = (args.eneEnd_-args.eneStart_) / QSfloat(args.eneSteps_)
    ene = args.eneStart_ + args.eneComplex_*1.0j
    sMats = {}
    for i in range(0,args.eneSteps_+1,1):
        if anaSmat is not None:
            anaSmat.setEnergy(ene)
            sMats[ene] = anaSmat.getMatrix()
        else:
            sMats[ene] = None
        ene += dEne
    return sMats 

def getRatSmat(args, anaSmat, anakCal, fitkCal, suppressCmdOut=False):
    sMats = getDiscreteAnaSmats(args, anaSmat)  
    return RatSMat(sMats, fitkCal, resultFileHandler=_getFileHandler(args, anakCal, fitkCal), suppressCmdOut=suppressCmdOut)

def getSmats(args, anakCal, fitkCal):
    anaSmat = getAnaSmat(args, anakCal)
    ratSmat = getRatSmat(args, anaSmat, anakCal, fitkCal)
    return (anaSmat, ratSmat)

def getDecimatedRatSmat(args, smats, anakCal, fitkCal, N, anaSmat=None, suppressCmdOut=False): #To get ratSmat for same data set as for getPolyRoots to allow comparison with Muller
    resultFileHandler = _getFileHandler(args, anakCal, fitkCal)
    decimator = Decimator(0, len(smats)-1, 0, resultFileHandler)
    newSMats, decStr = decimator.decimate(smats, N)
    if anaSmat is not None:
        popSmat(anaSmat, newSMats, smats)
    return RatSMat(newSMats, fitkCal, resultFileHandler=resultFileHandler, suppressCmdOut=suppressCmdOut)

def getPolyRoots(args, anakCal, fitkCal, resultsPath, mode, cmpValue=None):
    anaSmat = getAnaSmat(args, anakCal)
    smats = getDiscreteAnaSmats(args)
    PoleFinder(smats, fitkCal, _getFileHandler(args, anakCal, fitkCal), ENEFACTOR, 0, len(smats)-1, 0, args.distFactor_, args.cfSteps_, cmpValue=cmpValue, mode=mode, populateSmatCB=lambda sm1,sm2: popSmat(anaSmat, sm1, sm2), zeroValExp=args.zeroValExp_)

def popSmat(anaSmat, smats1, smats2=None):
    for ene in smats1:
        if smats1[ene] is None:
            anaSmat.setEnergy(ene)
            mat = anaSmat.getMatrix()
            smats1[ene] = mat
            if smats2 is not None:
                smats2[ene] = mat     
 
def dokSignIt(args, anaSignList, fitSignList, ratSignList, anaFun, ratFun, suppressCmdOut=False, signsAsList=False):
    try:
        first = True
        for anaSigns in anaSignList:
            if not first:
                print "\n"
            first = False
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