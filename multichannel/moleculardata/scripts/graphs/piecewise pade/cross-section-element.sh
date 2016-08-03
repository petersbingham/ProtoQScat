# This plots specified element of the cross-section matrix (-1 -1 in the parameters means a piecewise fit)

N=20
ROW=0
COL=0

cd ../../..
python ../qscat/numerical/run_plotcrosssection_element.py $N -1 -1 $ROW $COL
read -n1 -r -p "Press any key to continue..." key