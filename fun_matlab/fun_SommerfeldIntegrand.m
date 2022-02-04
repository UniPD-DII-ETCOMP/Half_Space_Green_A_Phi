%% brief Function to assemble the integrand for a given component of the MGF
function [y] = fun_SommerfeldIntegrand(rs,rt,e,k,freq,rho,krho,component,order)
mu0=4*pi*1e-7;
%% Compute spectral MGF
%krho=fun_SetRadialWaveNumber(krho); non serve a un cazzo
K=fun_HandleDyadicGreen(component,rs,rt,k,freq,e,mu0);
% kz1=mySqrtNew(k(1).^2-krho.^2);
% K=@(krho) K(krho)./kz1;
% Compute Bessel function of first kind
Jv=besselj(order,krho*rho);
%%
y=Jv.*K(krho).*(krho.^(order+1))/(2*pi);
%y=K(krho);
end