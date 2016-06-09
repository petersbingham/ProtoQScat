import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../../Utilities')
import general.simpplot as sp

xs = [3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0]
ys = [0.00100916651423,0.00033452148126,0.00064784098685,0.00063865709919,0.00064964557524,0.00065278686934,0.00065273417697,0.00065284102581,0.00065274916144,0.00065278530883,0.00065279744940,0.00065279742736,0.00065279742700]
sp.plotSingle("Imag values of Selected Roots Close to Resonance with df=0.01%, Emax=0.1194Ryds", xs, [ys], "N", "Imag Energy (Rydbergs)", path="Results/Values_ManPlot_Imag.png",markerSz=8,markWithLine=True,logy=False)