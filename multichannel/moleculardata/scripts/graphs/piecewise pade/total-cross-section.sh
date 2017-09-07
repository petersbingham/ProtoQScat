# This plots the total cross section

N=20

cd ../../..
python ../qscat/numerical/run_plottotalcrosssection.py $N
read -n1 -r -p "Press any key to continue..." key