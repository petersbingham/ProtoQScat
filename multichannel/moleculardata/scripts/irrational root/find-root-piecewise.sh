# This searches for a single root in a poecewise fit given a starting position (-1 -1 in the parameters means a piecewise fit).

SET_N=20
STARTING_POINT=1.0+0.0j

cd ../..
python ../qscat/numerical/run_findroot.py $SET_N -1 -1 $STARTING_POINT
read -n1 -r -p "Press any key to continue..." key