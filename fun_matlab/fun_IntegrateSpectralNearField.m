%% Integrate Near-Field region
%brief function to integrate the spectral MGF for values of krho in the
%spectral near-field, using a rectangular path deformation to avoid
%singularities
function [res] = fun_IntegrateSpectralNearField(rs,rt,e,k,freq,rho,component,order,a)
k_max=max(abs(k)); %wavenumbers
h=1e-2*abs(k_max);
%% Sommerefld Integrand
%% MGF integrand for a fixed imaginary part in the complex plane
%function handle
%f_real = @(krho) fun_IntegrandFixKrhoParts(real(krho),krho_fixed,'im',rs,rt,e,k,freq,rho,component,order);
%% MGF integrand for a fixed real part in the complex plane
%f_imag = @(krho) fun_IntegrandFixKrhoParts(imag(krho),krho_fixed,'re',rs,rt,e,k,freq,rho,component,order);
%%
tol=1e-6;
zero_offset=0;%1.0e-3; %0.0e-3 in strata????
%% First Path
fixed_part=zero_offset;
f_imag = @(x) fun_IntegrandFixKrhoParts(x,fixed_part,'re',rs,rt,e,k,freq,rho,component,order);
first_path=fun_GaussKronrodBoost(f_imag,0.0,h,tol);
%% Second Path
fixed_part=h;
f_real = @(x) fun_IntegrandFixKrhoParts(x,fixed_part,'im',rs,rt,e,k,freq,rho,component,order);
second_path=fun_GaussKronrodBoost(f_real,zero_offset,a,tol);
%% Third Path
fixed_part=a;
f_imag = @(x) fun_IntegrandFixKrhoParts(x,fixed_part,'re',rs,rt,e,k,freq,rho,component,order);
third_path=fun_GaussKronrodBoost(f_imag,h,0.0,tol);
%%
res=first_path+second_path+third_path;
end