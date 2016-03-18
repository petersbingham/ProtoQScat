from Base import *
from Table import *

def makeTable(df, index):
    filePath = base+"\\iMaxVar\\Results\\df="+str(df)+"_kmax="+str(index)+".txt"
    tabulateData(df, INDEX_EMAX, index, filePath)
    
makeTable(0.0001, 320)
makeTable(0.01, 192)
makeTable(0.01, 576)