# This plots the total cross section

N=20
SET_START=0
SET_END=1200
PLOT_START=0.001
PLOT_END=1.0
COMPLEXOFFSET=0.0
NUM_PLOT_POINTS=1500

cd ../../..
python ../qscat/numerical/run_ploteigenphase.py $N $SET_START $SET_END $PLOT_START $PLOT_END $COMPLEXOFFSET $NUM_PLOT_POINTS