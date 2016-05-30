STEP_TO_RYD = 0.000367408

def stepToRydStr(stepNum,pres=4):
    forStr = "{:."+str(pres)+"f}"
    return forStr.format(STEP_TO_RYD * float(stepNum))