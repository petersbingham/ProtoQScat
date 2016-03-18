import sys
import os

base_kvar =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base_kvar+'./../../..')

from Dat import *

p = Printer(True, False)
origstdout = sys.stdout

def tabulate(a,V,l,kmin,kmax,df):
    fileName = "/Results/a="+str(a)+" V="+str(V)+" l="+str(l)+" kmin="+str(kmin)+" kmax="+str(kmax)+" df="+str(df)+".txt"
    sys.stdout = origstdout
    d = getData(base_kvar+'/../../../Data/'+str(df)+'/'+str(a)+'/')
    sys.stdout = open(base_kvar+fileName, 'w')
    p.tabulatePoleData(d.getData(),a=a,V=V,l=l,kmin=kmin,kmax=kmax,df=df)
    sys.stdout.close()

tabulate(1.0,3.4,0,0.01,15.0,1.0)
  
for df in [0.001, 0.01, 0.1, 1.0]:
  for V in [1.3, 2.0, 3.4, 5.5, 7.6, 9.4, 9.7, 10.0, 10.3]:
    tabulate(1.0,V,0,0.01,15.0,df)

for df in [0.001, 0.01, 0.1, 1.0]: 
  for V in [0.4,1.0,3.0]:
    tabulate(2.0,V,0,0.01,3.5,df)