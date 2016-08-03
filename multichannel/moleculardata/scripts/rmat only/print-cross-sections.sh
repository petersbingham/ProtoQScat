# This calculates and writes to file values of the total cross section at the grid points (in eV). File at ./../../res.txt 

cd ../..
python ../qscat/numerical/run_printanalyticalcrosssections.py
read -n1 -r -p "Press any key to continue..." key