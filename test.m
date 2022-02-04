clear 
close all
clc
%%
eps0=8.8541878128e-12;
mu0=4*pi*1e-7;
c0=299792458;
format short
%%
%% BEGIN USER SETTINGS
%%
%% Frequency
freq=60.3e6;
%% Source and target point
l=c0/freq/30;
% source
rs=[+l -l/3 +l];
% target
rt=[-l/2 -l/4 -l];
%% Layer property
sigma=2; 
epsr=3;
%%
%% END USER SETTINGS
%%
epsr_eq=epsr-1j*sigma/(2*pi*freq*eps0);
%% Folders
dad=pwd; cd('fun_matlab'); addpath(genpath(pwd)); cd(dad)
dad=pwd; cd('fun_fortran'); addpath(genpath(pwd)); cd(dad)
%% Constants
lambda0=c0/freq;
omega=2*pi*freq;
k0=omega*sqrt(mu0*eps0);
eta0=sqrt(mu0/eps0);
%% Wavenumbers
k(1)=2*pi*freq*sqrt(eps0*mu0);
k(2)=2*pi*freq*sqrt(epsr_eq*eps0*mu0);
e=[eps0,eps0*epsr_eq];
%% DYADIC GREEN 
% matlab
[KA_m,Kphi_m]=fun_DyadicGreen(rs,rt,e,k,freq);
KA_m=KA_m*mu0;
Kphi_m=Kphi_m/eps0;
% fortran
try
    [KA_re,KA_i,Kphi_r,Kphi_i]=fun_DyadicGreen_f90(rs,rt,real(e),imag(e),real(k),imag(k),freq);
    KA_f=mu0*(KA_re+1j*KA_i);
    Kphi_f=(Kphi_r+1j*Kphi_i)/eps0;
catch
    KA_f=nan;
    Kphi_f=nan;
    warning('- MEX function not supported, try to re-mex it: run /fun_fortran/make.m.')
end
%%
disp('MATLAB')
KA_m
Kphi_m
disp('FORTRAN')
KA_f
Kphi_f
