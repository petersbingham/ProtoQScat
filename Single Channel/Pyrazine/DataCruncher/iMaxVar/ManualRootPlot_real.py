import sys
import os
base = os.path.dirname(os.path.realpath(__file__))
sys.path.insert(0,base+'/../../../../Utilities')
import General.SimpPlot as sp

xs = [3.0,4.0,5.0,6.0,7.0,8.0,9.0,10.0,11.0,12.0,13.0,14.0,15.0]
ys = [0.09127081964164,0.07780137160144,0.07714598042590,0.07680143868106,0.07693999245180,0.07693112841080,0.07693102556199,0.07693119716895,0.07693077955548,0.07693125892433,0.07693128907091,0.07693128907534,0.07693128908032]
sp.plotSingle("Real values of Selected Roots Close to Resonance with df=0.01%, Emax=0.1194Ryds", xs, [ys], "N", "Real Energy (Rydbergs)", path="Results/Values_ManPlot_Real.png",markerSz=8,markWithLine=True,logy=False)