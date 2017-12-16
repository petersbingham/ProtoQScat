# This returns all elastic roots by solving the polynomial.

N=20
SET_START=0
SET_END=1199

cd ../../..
python ../qscat/numerical/run_findpolyroots.py $N $SET_START $SET_END
#python -u ../qscat/numerical/run_findpolyroots.py $N $SET_START $SET_END > "./scripts/rational roots/roots/out.dat"
#python -m cProfile ../qscat/numerical/run_findpolyroots.py $N $SET_START $SET_END
#python -u -m cProfile ../qscat/numerical/run_findpolyroots.py $N $SET_START $SET_END > "./scripts/rational roots/roots/out.dat"
read -n1 -r -p "Press any key to continue..." key