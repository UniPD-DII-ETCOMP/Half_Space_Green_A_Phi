
# HALF SPACE GREEN FUNCTION
This directory contains a simplified matlab and fortran version of the STRATA library https://github.com/modelics/strata for the computation of the Green's function of an Half space background. It can be used for research and academic purposes. In this MATLAB/Fortran version only the half space case (two semi-infinite layers) is considered. Moreover, only the "integrate" technique is implemented and we only focus on the components of the dyadic K(A), and the scalar term K(φ), as defined in formulation-C of [1]
The original C++ Strata library, instead, also implements the multi-layer version of the Green's functions and other methods (DCIM, quasistatic) and features are provided. See https://github.com/modelics/strata for more details. 

-------------------------------------------------------------------

# References
[1] S. Sharma and P. Triverio, Strata: An Open-Source C++ Library for Computing Green's Functions for Layered Media, 2021 IEEE International Symposium on Antennas and Propagation and USNC-URSI Radio Science Meeting, 2021, Singapore
	
[2] K. A. Michalski and D. Zheng, "Electromagnetic scattering and radiation by surfaces of arbitrary shape in layered media. I. theory," IEEE Trans. Antennas Propag., vol. 38, no. 3, pp. 335–344, Mar. 1990.

[3] K. A. Michalski and J. R. Mosig, "Multilayered media Green's functions in integral equation formulations," IEEE Trans. Antennas Propag., vol. 45, no. 3, pp. 508–519, Mar. 1997.

[4] E. Simsek, Q. H. Liu, and B. Wei, "Singularity Subtraction for Evaluation of Green's Functions for Multilayer Media," IEEE Trans. Microw. Theory Tech., vol. 54, pp. 216–225, Jan 2006.

[5] K. A. Michalski and J. R. Mosig, "Efficient computation of Sommerfeld integral tails - methods and algorithms," J. Electromagn. Waves Appl., vol. 30, no. 3, pp. 281–317, 2016.

[6] M. I. Aksun, "A Robust Approach for the derivation of closed-form Green's functions," IEEE Trans. Microw. Theory Tech., vol. 44, pp. 651–658, May 1996.	

-------------------------------------------------------------------

# Note
The Matlab functions are not optimized. 
Small discrepancies between Matlab and Fortran functions can be obtained due to the numeric integrations. 

This toolbox contains parts or modifications of the following codes:

[1] STRATA, https://github.com/modelics/strata (MATLAB & FRTRAN porting)

[2] Jason Nicholson (2022). Bessel Zero Solver (https://www.mathworks.com/matlabcentral/fileexchange/48403-bessel-zero-solver), MATLAB Central File Exchange. Retrieved February 4, 2022.

[3] TOMS644, Bessel Functions of Complex Argument and Nonnegative Real Order. https://people.sc.fsu.edu/~jburkardt/f77_src/toms644/toms644.html

[4] quadpack, http://www.netlib.org/quadpack/

-------------------------------------------------------------------

# Contacts
Riccardo Torchio (riccardo.torchio@unipd.it)
Francesco Lucchini (francesco.lucchini@igi.cnr.it)
