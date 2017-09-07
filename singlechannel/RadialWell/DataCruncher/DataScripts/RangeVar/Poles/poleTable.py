import sys
import os
import argparse


parentArgs = argparse.ArgumentParser(description="Single Channel Radial Well. Data Analysis.")
parentArgs.add_argument("_a", help=".", type=float)
parentArgs.add_argument("_V", help=".", type=float)
parentArgs.add_argument("_l", help=".", type=int)
parentArgs.add_argument("_kmin", help=".", type=float)
parentArgs.add_argument("_kmax", help=".", type=float)
parentArgs.add_argument("_df", help=".", type=float)
args = parentArgs.parse_args()


base_kvar =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base_kvar+'./../../..')

from Dat import *

d = getData(base_kvar+'/../../../Data/'+str(args._df)+'/'+str(args._a)+'/')
p = Printer(True, False)

sys.stdout = open(base_kvar+"/Results/a="+str(args._a)+" V="+str(args._V)+" l="+str(args._l)+" kmin="+str(args._kmin)+" kmax="+str(args._kmax)+" df="+str(args._df)+".txt", 'w')
p.tabulatePoleData(d.getData(),a=args._a,V=args._V,l=args._l,kmin=args._kmin,kmax=args._kmax,df=args._df)