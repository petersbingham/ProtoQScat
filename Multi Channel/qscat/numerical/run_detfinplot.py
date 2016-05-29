from RatSMatWrap import *
import General.SimpPlot as sp
from runargsplot import *

kmat = RatSMatWrap("fort.19",args.subN_,args.subStart_,args.subEnd_)
vals = kmat.getFinDetRange(args.plotStart_, args.plotEnd_, args.plotComplex_, args.steps_)

sp.plotSingle("Pyrazine : det(Fin)", vals[0], [vals[1],vals[2]], 'Electron Energy (rydbergs)', 'det(Fin)', ['real','imag'])
