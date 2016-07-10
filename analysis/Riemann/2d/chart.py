import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../Utilities/General')
import simpplot as sp
sp.setSubPlotParameters(left=0.12, bottom=0.16, right=0.96, top=0.9, wspace=0.2, hspace=0.2)
sp.setImgSize(6,3)

P_REALS = [2.1,2.2,2.0,1.8,0.0,0.0,0.0]
P_IMAGS = [-0.2,-0.4,-0.6,-0.8,1.3,1.8,2.0]

E_REALS = [3.8,4.3,3.2,2.1,-1.5,-3.0,-4.0]
E_IMAGS = [-1.0,-2.0,-2.6,-3.1,0.0,0.0,0.0]

sp.plotSingle("p plane",P_REALS,[P_IMAGS],"Real (Ryd)","Imag",path="Results/pplane.png",markerSz=6,drawAxes=True)
sp.plotSingle("E plane",E_REALS,[E_IMAGS],"Real","Imag",path="Results/eplane.png",markerSz=6,drawAxes=True)