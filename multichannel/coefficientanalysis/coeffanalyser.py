from os import walk
import sys
import numpy as np

#Purpose of this is to compare two sets of coefficients by dividing one by the other (element wise). We can then easily see if one is -1 times the other (for eg).

INDEX_A = 0
INDEX_B = 1

mats = []
sets = []

def _getSectionStr(fitNum, coeff):
    return "\n@<*****************>@\nFit:" + str(fitNum) + "   Coeffs:" + str(coeff) + "\n@<*****************>@\n\n"

for setIndex in range(len(sys.argv)-1):
    lines = []
    sets = sys.argv[1:]
    setName = sets[setIndex]
    mats.append([])

    path = "./../CoefficientFiles/" + setName + "/"
    
    files = []
    for (dirpath, dirnames, filenames) in walk(path):
        files.extend(filenames)
        break
    files = sorted(files)
    
    lastSectionStr = ""
    
    for file in files:
        paras = file.split('.')[0].split('_')
        para_fit = int(paras[1])
        para_coeffNum = int(paras[2])
        para_coeffType = INDEX_A if paras[0]=="A" else INDEX_B
        
        while len(mats[setIndex]) <= para_fit:
            mats[setIndex].append([])
        if len(mats[setIndex][para_fit]) == 0:
            mats[setIndex][para_fit].append([])
            mats[setIndex][para_fit].append([])
        mats[setIndex][para_fit][para_coeffType].append(np.loadtxt(path+file, dtype=np.complex128, delimiter=","))
        
        lines.append("\n")    
        sectionStr = _getSectionStr(paras[1], paras[0])
        if sectionStr != lastSectionStr:
            lines.append(sectionStr)
            lastSectionStr = sectionStr
    
        for line in open(path+file, 'r'):    
            lines.append(line)
        
    file = open("./Results/Coeffs_" + str(setIndex) + "_" + setName+".dat", 'w')
    for line in lines:
        file.write(line)
    file.close()

if len(mats) > 1:
    for setIndex1 in range(0,len(mats)):
        for setIndex2 in range(setIndex1+1,len(mats)):
            file = open("./Results/CoeffCmp_" + str(setIndex1) + "_" + str(setIndex2) + "_" + setName+".dat", 'w')
            for fitIndex in range(len(mats[setIndex1])):
                for typeIndex in range(len(mats[setIndex1][fitIndex])):
                    sectionStr = _getSectionStr(fitIndex, "A" if typeIndex == INDEX_A else "B")
                    file.write(sectionStr)
                    for matIndex in range(len(mats[setIndex1][fitIndex][typeIndex])):
                        coeffMat = mats[setIndex1][fitIndex][typeIndex][matIndex]
                        cmpCoeffMat = mats[setIndex2][fitIndex][typeIndex][matIndex]
                        with np.errstate(divide='ignore', invalid='ignore'):
                            file.write(str(coeffMat/cmpCoeffMat)+"\n\n")
    