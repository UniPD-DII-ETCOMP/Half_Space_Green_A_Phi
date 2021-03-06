clear 
close all
clc

mex -O -v -largeArrayDims -output fun_DyadicGreen_f90 ...
    myfargs.f90 ...
    mySqrtNew.f90 ...
    fun_BesselJ_NextZero.f90 ...
    fun_IntegrateSpectralFarField.f90 ...
    fun_IntegrateSpectralNearField.f90 ...
    fun_LevinSidi.f90 ...
    fun_PartExtrap.f90 ...
    fun_PartSum.f90 ...
    fun_TanhSinh.f90 ...
    fun_Term.f90 ...
    fun_TruncIndex.f90 ...
    fun_GA1_ii_Im.f90 ...
    fun_GA1_ii_Re.f90 ...
    fun_GA1_mi_Im.f90 ...
    fun_GA1_mi_Re.f90 ...
    fun_GA2_ii_Im.f90 ...
    fun_GA2_ii_Re.f90 ...
    fun_GA2_mi_Im.f90 ...
    fun_GA2_mi_Re.f90 ...
    fun_GA3_ii_Im.f90 ...
    fun_GA3_ii_Re.f90 ...
    fun_GA3_mi_Im.f90 ...
    fun_GA3_mi_Re.f90 ...
    fun_GA4_ii_Im.f90 ...
    fun_GA4_ii_Re.f90 ...
    fun_GA4_mi_Im.f90 ...
    fun_GA4_mi_Re.f90 ...    
    fun_Gphi_ii_Im.f90 ...
    fun_Gphi_ii_Re.f90 ...
    fun_Gphi_mi_Im.f90 ...
    fun_Gphi_mi_Re.f90 ...    
    fun_GaussKronrodBoost.f90 ...
    fun_ComputeMGF_Integration_mex.f90 fun_ComputeMGF_Integration.f90 ...
    mytoms644.f ...
    dqags.f d1mach.f dqagse.f dqelg.f dqk21.f ...
    dqpsrt.f fdump.f i1mach.f j4save.f xercnt.f ...
    xerhlt.f xermsg.f xerprn.f xersve.f xgetua.f