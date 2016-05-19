ProtoQScat
==========

Protoype for quantum scattering library. Main two folders are Utilities and Multi Channel. Single Channel folder contains some data analysis code used with results from another project. The Analysis folder contains an examination of the Riemann energy surfaces.

<br />
Dependencies
____________
 * python 2.7
 * numpy
 * scipy
 * sympy
 * matplotlib

The main routines can be found in Multi Channel\RatSMat.py. The pole location routines are in Multi Channel\PoleFinder.py.
 
<br /> 
Result Archiving
________________

Depending on the calculation configuration significant time may be required to complete the calculations. For this reason and to assist with the organisation of results, generated files are automatically stored according to the structure below.

 * Dir: Fit Type and Set ID String.
   * Dir: Method and parameters used to calculated the coefficients.
     * Dir: Nature of the fit (single or piecewise).
       * Dir: Coefficient Files.
         * Dir: Number of Points, start and end indices.
           * Files: One for each coefficient matrix
         * Dir: Nature of the method and parameters used to calculate the roots.
           * Dir: Roots
             * Files: One for each number of Points, start and end indices.
           * Dir: Poles, according to technique (doubling or incrementing N) and the threshold comparison factor.
             * Files: One for each number of Points, start and end indices.
  
The software will attempt to load any prior dependencies before doing calculations. So, for example, if we have already calculated a set of poles and then change only the threholds, the software will load the previously calculated roots as a starting point, rather than recalculating the coefficients and roots again. Similarly, if we've already calculated a set of roots and then only change the method or parameters for the root calculation the prior coefficients will be used as a starting point. 

These files are by default, stored in the Results Folder in the "Multi Channel" directory. 

<br />
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

<br />
Pyrazine Tests
______________

Based on scattering data obtained from Pyrazine molecule. Scripts can be found in the following folders:
 * Multi Channel\Pyrazine\Analytical
 * Multi Channel\Pyrazine\Comparison
 * Multi Channel\Pyrazine\PieceWise
 * Multi Channel\Pyrazine\PoleFinders
 * Multi Channel\Pyrazine\Single
 
