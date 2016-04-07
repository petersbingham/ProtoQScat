ProtoQScat
==========

Protoype for quantum scattering library. Main two folders are Utilities and Multi Channel. Single Channel folder contains some data analysis code used with results from another project. The Analysis folder contains an examination of the Riemann energy surfaces.

Dependencies
____________
 * python 2.7
 * numpy
 * scipy
 * sympy
 * matplotlib

The main routines can be found in Multi Channel\RatSMat.py. The pole location routines are in Multi Channel\PoleFinder.py.
 
Multi-channel Analytical Tests
______________________________

These are based on the analytical solutions described in: R G Newton, Scattering Theory of Waves and Particles, Springer-Verlag, New York, 2nd edition (1982).
Results scripts (cross sections, find roots etc) can be found in the following folders. Simply execute them:
 * Multi Channel\TwoChannel\Bats_1_2_2_0_0_0
 * Multi Channel\TwoChannel\Bats_1_2_2_0_0_1
 * Multi Channel\TwoChannel\Bats_1_2_2_0_2_0
 * Multi Channel\TwoChannel\Bats_1_2_2_0_2_1  

The numbers refer to the parameters of the well: _width_depth1_depth2_threshold1_threshold2_couplingfactor. The numbers in the script names refer to the fit parameters.

Further scripts can be found in:
 * TwoChannel\Analytical
These will obtain results purely from the analytical solutions.


Pyrazine Tests
______________

Based on scattering data obtained from Pyrazine molecule. Scripts can be found in the following folders:
 * Multi Channel\Pyrazine\Analytical
 * Multi Channel\Pyrazine\Comparison
 * Multi Channel\Pyrazine\PieceWise
 * Multi Channel\Pyrazine\PoleFinders
 * Multi Channel\Pyrazine\Single
 
