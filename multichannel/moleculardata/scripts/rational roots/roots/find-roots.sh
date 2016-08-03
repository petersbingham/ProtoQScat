# This returns all elastic roots by solving the polynomial.

N=20
SET_START=0
SET_END=2000

cd ../../..
python ../qscat/numerical/run_findpolyroots.py $N $SET_START $SET_END
read -n1 -r -p "Press any key to continue..." key