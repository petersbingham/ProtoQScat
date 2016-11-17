# Locates poles by comparing roots for N=N*2. Positive signs are used throughout.

SET_START=0
SET_END=-1 # To end of data
SET_OFFSET=0

#delimit CF_STEPS and DIST_FACTOR with commas for several values to loop over
CF_STEPS=1,2
DIST_FACTOR=0.1,0.01,0.001,0.0001,0.00001,0.000001,0.0000001,0.00000001,0.000000001,0.0000000001,0.00000000001,0.000000000001,0.0000000000001
ZERO_VALUE_EXPONENT=7

N_MIN=4
N_MAX=40

FLAG_RMATRIX_POLES_INDEX=0

cd ../../..
python ../qscat/numerical/polecalculatorwrap.py $SET_START $SET_END $SET_OFFSET 4 $CF_STEPS $DIST_FACTOR $ZERO_VALUE_EXPONENT $N_MIN $N_MAX $FLAG_RMATRIX_POLES_INDEX
read -n1 -r -p "Press any key to continue..." key