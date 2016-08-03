# This searches for a single root in a single rational fit from a given starting guess and its conjugate.

N=20
SET_START=0
SET_END=2000
STARTING_POINT=1.0+0.0j

cd ../..
python ../qscat/numerical/run_findroot_withconj.py $N $SET_START $SET_END $STARTING_POINT
read -n1 -r -p "Press any key to continue..." key