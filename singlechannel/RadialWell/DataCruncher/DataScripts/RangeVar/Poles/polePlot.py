import sys
import os
import argparse

'''
parentArgs = argparse.ArgumentParser(description="Single Channel Radial Well. Data Analysis.")
parentArgs.add_argument("a_", help="Width of Well.", type=float)
parentArgs.add_argument("results_", help="1-Root Table, 2-Pole Table, 3-Root Total Plot, 4-Pole Total Plot, 5-Root Scat Plot, 6-Pole Scat Plot", nargs='?', type=int, default=1)
parentArgs.add_argument("kmin_", help="Value of kmin", nargs='?', type=float, default=0.01)
parentArgs.add_argument("chopEnd_", help="For scattering number of spurious states to remove.", nargs='?', type=int, default=None)
parentArgs.add_argument("logx_", help="Plot logarithmic x axis.", nargs='?', type=bool, default=False)
args = parentArgs.parse_args()
'''

a = 2.0
df = 0.001
base_kvar =  os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base_kvar+'/../../..')

from Dat import *

import general.simpplot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.12, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(5,4)

d = getData(base_kvar+'/../../../Data/'+str(df)+'/'+str(a)+'/')
p = Printer(True, False)
xs, ys = p.plotPoles_E(d.getData(),a=a,V=1.0,l=0,kmin=0.01,kmax=4.0,df=df)

sp.plotScat("", xs, ys)
      
    
  
  
  
  
  
