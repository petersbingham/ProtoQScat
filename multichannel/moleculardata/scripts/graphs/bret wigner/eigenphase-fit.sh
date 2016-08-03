# This fits to poles as specified in the sysdesc.py

ENE_START=0.0
ENE_END=1.0
NUM_PLOT_POINTS=5000
POLY_FIT_ORDER=2

cd ../../..
python ../qscat/numerical/run_ploteigenphasefit.py $ENE_START $ENE_END $NUM_PLOT_POINTS $POLY_FIT_ORDER
read -n1 -r -p "Press any key to continue..." key