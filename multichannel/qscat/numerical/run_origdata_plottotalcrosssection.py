from ratsmatwrap import *

try:
    kmat = RatSMatWrap(FILENAME)
    rmatxs = kmat.getTotalDiscreteXS()
   
    rmatxs.setDetails("R matrix", ['green'])
    
    sm.plotAll("Log", [rmatxs], title="Fitted Numerical Data Cross Sections", xlabel="Electron Energy (rydbergs)", ylabel="Cross Section (bohr^2)") 
except (sm.MatException) as inst:
    print str(inst)
    sys.exit()