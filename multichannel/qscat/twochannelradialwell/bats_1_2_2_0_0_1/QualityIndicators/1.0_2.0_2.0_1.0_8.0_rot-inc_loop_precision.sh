# To simulate results that had been obtained for the pyrazine molecular. The selected depth results in a resonance of similar width.

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

for ANALYTIC_TYPE in $(seq 2 2 18)
do
  python ./run_polecalculatorwrap.py 1.0 2.0 2.0 0.0 0.0 1.0 1.0 8.0 0.0 1000 $SET_START $SET_END $SET_OFFSET 3 $CF_STEPS $START_DIST_THRES $AMALG_THRES $ZERO_VALUE_EXPONENT $N_MIN $N_MAX $ANALYTIC_TYPE
done

read -n1 -r -p "Press any key to continue..." key