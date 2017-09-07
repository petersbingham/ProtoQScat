# This plots specified row of the cross-section matrix (-1 -1 in the parameters means a piecewise fit)

N=20
ROW=0

cd ../../..
python ../qscat/numerical/run_plotcrosssection_row.py $N -1 -1 $ROW
read -n1 -r -p "Press any key to continue..." key