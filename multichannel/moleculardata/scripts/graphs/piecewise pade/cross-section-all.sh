# This plots all elements of the cross-section matrix

N=20

cd ../../..
python ../qscat/numerical/run_plotcrosssection_all.py $N
read -n1 -r -p "Press any key to continue..." key