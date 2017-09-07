ProtoQScat
==========

Prototype for quantum scattering library. Main two folders are utilities and multichannel. singlechannel folder contains some data analysis code used with results from another project. The analysis folder contains an examination of the Riemann energy surfaces and various other investigations.

<br />
Dependencies
____________
 * python 2.7
 * numpy+mkl
 * scipy
 * sympy
 * matplotlib
 * tabulate

The main routines can be found in multichannel/qscat/ratsmat/__init__.py. The pole location routines are in multichannel/qscat/ratsmat/polefinder.py.
 
<br />
Multi-channel Analytical Tests
______________________________

These are based on the analytical solutions described in: R G Newton, Scattering Theory of Waves and Particles, Springer-Verlag, New York, 2nd edition (1982).
Results scripts (cross sections, find roots etc) can be found in the following folders. Simply execute them:
 * multichannel/qscat/twochannelradialwell/bats_1_2_2_0_0_0
 * multichannel/qscat/twochannelradialwell/bats_1_2_2_0_0_1
 * multichannel/qscat/twochannelradialwell/bats_1_2_2_0_2_0
 * multichannel/qscat/twochannelradialwell/bats_1_2_2_0_2_1  

The numbers refer to the parameters of the well: _width_depth1_depth2_threshold1_threshold2_couplingfactor. The numbers in the script names refer to the fit parameters.

Further scripts can be found in:
 * multichannel/qscat/twochannelradialwell/analytical
These will obtain results purely from the analytical solutions.

<br />
Molecular Data
______________

To add new data for analysis, follow file format as illustrated in existing examples (eg multichannel/moleculardata/rmatrixdata/pyrazine_3ch.19) or create your own version of multichannel/moleculardata/rfortmatreader.py:
 * In multichannel/moleculardata/rmatrixdata add the data file.
 * In the multichannel/moleculardata/sysdesc.py file add a code block and the details describing your data. Follow the convention in the existing examples.
 * Various scripts can be found here: multichannel/moleculardata/scripts. Open the script to see documentation and varible parameters.
 * Additional lengthy parameters can be adjusted in the file: scriptparameters.py.
 
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