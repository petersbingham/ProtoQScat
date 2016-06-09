from ratsmatwrap import *
import general.simpplot as sp
from runargsplot import *

kmat = RatSMatWrap(FILENAME,args.subN_,args.subStart_,args.subEnd_)
vals = kmat.getFinDetRange(args.plotStart_, args.plotEnd_, args.plotComplex_, args.steps_)

sp.plotSingle("Numercal Data : det(Fin)", vals[0], [vals[1],vals[2]], 'Electron Energy (rydbergs)', 'det(Fin)', ['real','imag'])
