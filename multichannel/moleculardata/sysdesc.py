DESC_STR = "Pyrazine Elastic 3ch"

"Urazil App Elastic Equil SE 6ch"

"Urazil Ap Elastic Equil SE 10ch"
"Urazil Ap Elastic Equil SE 10ch v2"
"Urazil Ap Elastic 0.1 SE 10ch"
"Urazil Ap Elastic 0.2 SE 10ch"

"Urazil Ap Elastic Equil SEP 10ch"

"para-benzoquinone B2u 3ch"
"para-benzoquinone se 3ch"

"CO2 a1 9ch"

"HCOOH se 15ch"

if DESC_STR == "Pyrazine Elastic 3ch":
    ARCHIVE_BASE_STR = "Pyrazine"
    FILENAME = "rmatrixdata/pyrazine_3ch.19"

    THRESHOLDS = [0.0]*3
    LS = [3,5,5]

elif DESC_STR == "Urazil App Elastic Equil SE 6ch":
    ARCHIVE_BASE_STR = "Urazil_6chanEleastic"
    FILENAME = "rmatrixdata/uracil_app_el_eq_se_6ch.19"
    
    THRESHOLDS = [0.0]*6
    LS = [1,2,2,3,3,3]
    
elif DESC_STR == "Urazil Ap Elastic Equil SE 10ch":
    ARCHIVE_BASE_STR = "Urazil_10chanEleastic"
    FILENAME = "rmatrixdata/uracil_ap_el_eq_se_10ch.19"
    
    THRESHOLDS = [0.0]*10
    LS = [0,1,1,2,2,2,3,3,3,3]

elif DESC_STR == "Urazil Ap Elastic Equil SE 10ch v2":
    ARCHIVE_BASE_STR = "Urazil_10chanEleastic_v2"
    FILENAME = "rmatrixdata/uracil_ap_el_eq_se_10ch_v2.19"
    
    THRESHOLDS = [0.0]*10
    LS = [0,1,1,2,2,2,3,3,3,3]
    
elif DESC_STR == "Urazil Ap Elastic 0.1 SE 10ch":
    ARCHIVE_BASE_STR = "Urazil_10chanEleastic_0.1"
    FILENAME = "rmatrixdata/uracil_ap_el_0.1_se_10ch.19"
    
    THRESHOLDS = [0.0]*10
    LS = [0,1,1,2,2,2,3,3,3,3]
    
elif DESC_STR == "Urazil Ap Elastic 0.2 SE 10ch":
    ARCHIVE_BASE_STR = "Urazil_10chanEleastic_0.2"
    FILENAME = "rmatrixdata/uracil_ap_el_0.2_se_10ch.19"
    
    THRESHOLDS = [0.0]*10
    LS = [0,1,1,2,2,2,3,3,3,3]
    
elif DESC_STR == "Urazil Ap Elastic Equil SEP 10ch":
    ARCHIVE_BASE_STR = "Urazil_10chanEleastic_SEP"
    FILENAME = "rmatrixdata/uracil_ap_el_eq_sep_10ch.19"
    
    THRESHOLDS = [0.0]*10
    LS = [0,1,1,2,2,2,3,3,3,3]
    
elif DESC_STR == "para-benzoquinone B2u 3ch":
    ARCHIVE_BASE_STR = "Para-benzoquinone_B2u_3ch"
    FILENAME = "rmatrixdata/para-benzoquinone_B2u_3ch.19"
    
    THRESHOLDS = [0.0]*3
    LS = [1,3,3]
    
elif DESC_STR == "para-benzoquinone se 3ch":
    ARCHIVE_BASE_STR = "para-benzoquinone_se_3ch"
    FILENAME = "rmatrixdata/para-benzoquinone_se_3ch.19"
    
    THRESHOLDS = [0.0]*3
    LS = [2,4,4]
    
elif DESC_STR == "CO2 a1 9ch":
    ARCHIVE_BASE_STR = "co2 a1"
    FILENAME = "rmatrixdata/co2_a1.19"
    
    THRESHOLDS = [0.0]*9
    LS = [0,1,2,2,3,3,4,4,4]
    
elif DESC_STR == "HCOOH se 15ch":
    ARCHIVE_BASE_STR = "hcooh se"
    FILENAME = "rmatrixdata/hcooh_se.19"
    
    THRESHOLDS = [0.0]*15
    LS = [0,1,1,2,2,2,3,3,3,3,4,4,4,4,4]
    
NUMCHANNELS = len(THRESHOLDS)
