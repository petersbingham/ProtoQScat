from sysdesc import *

# RMATRIX_POLES - Any pole positions as determined from the R-matrix calculations.
# CALCULATED_POLES - Any pole positions as determined from pade fit and analysis.

if DESC_STR == "Pyrazine Elastic 3ch":
    RMATRIX_POLES = [0.076641273-0.0006538903j]
    CALCULATED_POLES = [0.0766413036828500-0.0006529903921700j]

elif DESC_STR == "Urazil App Elastic Equil SE 6ch":
    RMATRIX_POLES = []
    CALCULATED_POLES = [0.1727679402440704-0.0089977530964865j, 
                        0.3408411506411183-0.0099343023220512j, 
                        0.6090283232752558-0.0759198056418777j, 
                        0.4118190571242367-0.2389125166701332j]
    
elif DESC_STR == "Urazil Ap Elastic Equil SE 10ch":
    RMATRIX_POLES = []
    CALCULATED_POLES = [0.6351845090672091-0.0113185404255686j, 
                        0.6827775357511159-0.0420373163357024j, 
                        0.8156992150572238-0.0477551269354513j, 
                        0.7895482191079259-0.0807704489650296j]

elif DESC_STR == "Urazil Ap Elastic Equil SE 10ch v2":
    RMATRIX_POLES = []
    CALCULATED_POLES = []
    
elif DESC_STR == "Urazil Ap Elastic 0.1 SE 10ch":
    RMATRIX_POLES = []
    CALCULATED_POLES = []
    
elif DESC_STR == "Urazil Ap Elastic 0.2 SE 10ch":
    RMATRIX_POLES = []
    CALCULATED_POLES = []
    
elif DESC_STR == "Urazil Ap Elastic Equil SEP 10ch":
    RMATRIX_POLES = []
    CALCULATED_POLES = []