from ratsmatwrap import *
import general.simpplot as sp
from runargsplot import *

kmat = RatSMatWrap(FILENAME,args.subN_,args.subStart_,args.subEnd_)
vals = kmat.getFinDetRealRange(args.plotStart_, args.plotEnd_, args.offset_, args.steps_)
dvals = kmat.getDiffFinDetRealRange(args.plotStart_, args.plotEnd_, args.offset_, args.steps_)

sp.plotSingle("Numercal Data : real det(Fin). Offset: " + str(args.offset_), vals[0], [vals[1],dvals[1]], 'Electron Energy (rydbergs)', 'det(Fin)', ['function','differential'])
sp.plotSingle("Numercal Data : imag det(Fin). Offset: " + str(args.offset_), vals[0], [vals[2],dvals[2]], 'Electron Energy (rydbergs)', 'det(Fin)', ['function','differential'])
