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
N_MAX=30

ANALYTIC_TYPE=0 #For the analytical calculation. 0 indicates to use the globally selected type (in qstype.py), -1 to simulate a 32-bit float, -2 for 64-bit. Positive number indicates the number of points to truncate the mantissa to.

function run {
  START="$(echo print $2-3.5 | python)"
  END="$(echo print $2+3.5 | python)"
  python ../../run_polecalculatorwrap.py 1.0 $1 $1 0.0 2.0 1.0 $START $END 0.0 1000 $SET_START $SET_END $SET_OFFSET 3 $CF_STEPS $START_DIST_THRES $AMALG_THRES $ZERO_VALUE_EXPONENT $N_MIN $N_MAX $ANALYTIC_TYPE
}

run 0 5.839
run 1 6.482
run 2 6.131
run 3 5.52
run 4 4.791
run 5 3.996
run 6 3.159
run 7 2.293
run 8 1.405
run 9 0.503
run 10 -0.412

: <<'END'
run 0 5.839
run 0.5 6.415
run 1 6.482
run 1.5 6.356
run 2 6.131
run 2.5 5.846
run 3 5.52
run 3.5 5.166
run 4 4.791
run 4.5 4.4
run 5 3.996
run 5.5 3.582
run 6 3.159
run 6.5 2.729
run 7 2.293
run 7.5 1.851
run 8 1.405
run 8.5 0.956
run 9 0.503
run 9.5 0.046
run 10 -0.412
END

read -n1 -r -p "Press any key to continue..." key