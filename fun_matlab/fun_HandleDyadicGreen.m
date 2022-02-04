function [y] = fun_HandleDyadicGreen(component,rs,rt,k,freq,e,mu0)
if component == 5
    %Scalar Potential Green function
    y = @(krho) fun_Kphi_spectral(krho,rs,rt,k,freq,e,mu0);
else
    %Tensorial A vector Green function
    y = @(krho) fun_KA_spectral(krho,rs,rt,k,freq,e,mu0,component);
end