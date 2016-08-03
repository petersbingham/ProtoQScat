# This plots the determinant of the Fin jost (-1 -1 in the parameters means a piecewise fit)

N=20
PLOT_START=0.02
PLOT_END=0.12
COMPLEXOFFSET=0.0
NUM_PLOT_POINTS=500

cd ../../..
python ../qscat/numerical/run_detfinplot.py $N -1 -1 $PLOT_START $PLOT_END $COMPLEXOFFSET $NUM_PLOT_POINTS
read -n1 -r -p "Press any key to continue..." key