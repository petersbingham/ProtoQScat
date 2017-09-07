# Locates poles by comparing roots for N=N*2. Routine runs for all sign permutations.

SET_START=0
SET_END=-1 # To end of data
SET_OFFSET=0

#delimit CF_STEPS and DIST_THRES with commas for several values to loop over
CF_STEPS=2,3
START_DIST_THRES=0.0001
AMALG_THRES=0.0001
ZERO_VALUE_EXPONENT=7

N_MIN=4
N_MAX=40

FLAG_RMATRIX_POLES_INDEX=0

cd ../../..
python ../qscat/numerical/polecalculatorwrap.py $SET_START $SET_END $SET_OFFSET 0 $CF_STEPS $START_DIST_THRES $AMALG_THRES $ZERO_VALUE_EXPONENT $N_MIN $N_MAX $FLAG_RMATRIX_POLES_INDEX
read -n1 -r -p "Press any key to continue..." key