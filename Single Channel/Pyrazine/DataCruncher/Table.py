from Dat import *

def tabulateData(df, pathIndex, pathValue, filePath):
    d = getData("./Data_"+str(df)+"/",pathIndex)
    p = Printer()
    sys.stdout = open(filePath, 'w')
    p.tabulatePoleData(d.getData(), pathValue)