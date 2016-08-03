# This searches for single roots for all Ns between 4 and 32 from a given starting guess.

SET_START=0
SET_END=1000
STARTING_POINT=1.0+0.0j

cd ../..
python ../qscat/numerical/run_findroot_inc_ns.py $SET_START $SET_END $STARTING_POINT
read -n1 -r -p "Press any key to continue..." key