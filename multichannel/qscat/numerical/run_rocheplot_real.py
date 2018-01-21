from ratsmatwrap import *
import general.simpplot as sp
from runargsplot import *

kmat = RatSMatWrap(FILENAME,args.subN_,args.subStart_,args.subEnd_)
vals = kmat.getRocheRealRange(args.plotStart_, args.plotEnd_, args.offset_, args.steps_)

sp.plotSingle("Numercal Data : roche for real energy. Offset: " + str(args.offset_), vals[0], [vals[1],vals[2]], 'Electron Energy (rydbergs)', 'roche', ['real','imag'])
