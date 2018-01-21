# This plots the determinant of the Fin jost

N=20
SET_START=0
SET_END=1000
PLOT_START=0.02
PLOT_END=0.12
COMPLEXOFFSET=0.0
NUM_PLOT_POINTS=500

cd ../../..
python ../qscat/numerical/run_detfinplot_imag.py $N $SET_START $SET_END $PLOT_START $PLOT_END $OFFSET $NUM_PLOT_POINTS
read -n1 -r -p "Press any key to continue..." key