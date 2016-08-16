ROOT_PATH="E:\\Peter's Documents\\PhD\\Code\\Git\\ProtoQScat\\multichannel\\qscat\\ratsmat\\Results"
SYS_PATH="Urazil_6chanEleastic_[1.0, 1.0, 1.0, 1.0, 1.0, 1.0]_0_1199"
COEFFS="COEFFS-mpmath_qr_solve_dps(norm None) DPS100"
FIT_TYPE="SingleFit"
ROOTS="ROOTS-v2_sympy_det(method berkowitz),sympy_Poly(),sympy_nroots(cleanup True, maxsteps 5000, n 100)"
POLES="Poles_incN_cfStep2_dk0.01_zk1e-07"

python ../qscat/ratsmat/ResultsAnalyser.py "$ROOT_PATH"\\"$SYS_PATH"\\"$COEFFS"\\"$FIT_TYPE" "$ROOTS" "$POLES" 
read -n1 -r -p "Press any key to continue..." key