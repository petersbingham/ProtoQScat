# Locates poles by comparing roots for N=N*2. Positive signs are used throughout.

SET_START=0
SET_END=-1 # To end of data
SET_OFFSET=0

#delimit CF_STEPS and DIST_THRES with commas for several values to loop over
CF_STEPS=2,3
START_DIST_THRES=0.0001
AMALG_THRES=0.0001
ZERO_VALUE_EXPONENT=7

N_MIN=4
N_MAX=32

ANALYTIC_TYPE=0 #For the analytical calculation. 0 indicates to use the globally selected type (in qstype.py), -1 to simulate a 32-bit float, -2 for 64-bit. Positive number indicates the number of points to truncate the mantissa to.

for DEPTH in $(seq 10.0 0.01 10.1)
do
  python ../../run_polecalculatorwrap.py 1.0 $DEPTH $DEPTH 0.0 0.0 1.0 1.0 8.0 0.0 1000 $SET_START $SET_END $SET_OFFSET 3 $CF_STEPS $START_DIST_THRES $AMALG_THRES $ZERO_VALUE_EXPONENT $N_MIN $N_MAX $RESULTS_TYPE $ANALYTIC_TYPE
done

read -n1 -r -p "Press any key to continue..." key