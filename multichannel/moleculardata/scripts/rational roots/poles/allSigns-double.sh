# Locates poles by comparing roots for N=N*2. Routine runs for all sign permutations.

SET_START=0
SET_END=1000
SET_OFFSET=0

CF_STEPS=1
DIST_FACTOR=0.1
ZERO_VALUE_EXPONENT=7

FLAG_RMATRIX_POLES_INDEX=0

cd ../../..
python ../qscat/numerical/polefinderwrap.py $SET_START $SET_END $SET_OFFSET 0 $CF_STEPS $DIST_FACTOR $ZERO_VALUE_EXPONENT $FLAG_RMATRIX_POLES_INDEX
read -n1 -r -p "Press any key to continue..." key