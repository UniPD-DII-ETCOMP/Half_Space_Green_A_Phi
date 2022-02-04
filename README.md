
# HALF SPACE GREEN FUNCTION

This directory contains a simplified matlab and fortran version of the STRATA library https://github.com/modelics/strata for the computation of the Green's function of an Half space background. 
It can be used for research and academic purposes. 
In this MATLAB/Fortran version only the half space case (two semi-infinite layers) is considered. Moreover, only the "integrate" technique is implemented and we only focus on the components of the dyadic K(A), and the scalar term K(φ), as defined in formulation-C of [1]
The original C++ Strata library, instead, also implements the multi-layer version of the Green's functions and other methods (DCIM, quasistatic) and features are provided. See https://github.com/modelics/strata for more details. 
-------------------------------------------------------------------

# References

@INPROCEEDINGS{strata,
	author={S. {Sharma} and P. {Triverio}},
	booktitle={2021 {IEEE} International Symposium on Antennas and Propagation and {USNC-URSI} Radio Science Meeting},
	title={Strata: An Open-Source {C++} Library for Computing {Green's} Functions for Layered Media},
	year={2021},
	month={Dec.},
	address = {Singapore}}
	
[1] K. A. Michalski and D. Zheng, "Electromagnetic scattering and radiation by surfaces of arbitrary shape in layered media. I. theory," IEEE Trans. Antennas Propag., vol. 38, no. 3, pp. 335–344, Mar. 1990.

[2] K. A. Michalski and J. R. Mosig, "Multilayered media Green's functions in integral equation formulations," IEEE Trans. Antennas Propag., vol. 45, no. 3, pp. 508–519, Mar. 1997.

[3] E. Simsek, Q. H. Liu, and B. Wei, "Singularity Subtraction for Evaluation of Green's Functions for Multilayer Media," IEEE Trans. Microw. Theory Tech., vol. 54, pp. 216–225, Jan 2006.

[4] K. A. Michalski and J. R. Mosig, "Efficient computation of Sommerfeld integral tails - methods and algorithms," J. Electromagn. Waves Appl., vol. 30, no. 3, pp. 281–317, 2016.

[5] M. I. Aksun, "A Robust Approach for the derivation of closed-form Green's functions," IEEE Trans. Microw. Theory Tech., vol. 44, pp. 651–658, May 1996.	


-------------------------------------------------------------------

# Note
 
The Matlab functions are not optimized. 
Small discrepancies between Matlab and Fortran functions can be obtained due to the numeric integrations. 

Contacts
-----------------------
Riccardo Torchio (riccardo.torchio@unipd.it)
Francesco Lucchini (francesco.lucchini@igi.cnr.it)