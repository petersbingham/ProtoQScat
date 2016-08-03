# This plots the total cross section

N=20
SET_START=0
SET_END=1000
PLOT_START=0.001
PLOT_END=0.6
COMPLEXOFFSET=0.0
NUM_PLOT_POINTS=500

cd ../../..
python ../qscat/numerical/run_plottotalcrosssection.py $N $SET_START $SET_END $PLOT_START $PLOT_END $COMPLEXOFFSET $NUM_PLOT_POINTS