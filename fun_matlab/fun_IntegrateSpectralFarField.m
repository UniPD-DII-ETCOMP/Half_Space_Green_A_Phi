%% Integrate Far-Field region
function [y] = fun_IntegrateSpectralFarField(rs,rt,e,k,freq,rho,component,order,a)
%% Initialization
tol=1e-6;
b=fun_BesselJ_NextZero(order,a*rho)/rho;
q=pi/rho;
%% MGF Integrand function
f=@(krho) fun_SommerfeldIntegrand(rs,rt,e,k,freq,rho,krho,component,order);
%% For small rho, Bessel function oscilllations are slow, so use quadrature
if rho < 1e-15
    result=fun_GaussKronrodBoost(f,a,Inf,tol);
else
    result=fun_PartExtrap(f,b,q,tol);
end
bridge=0.0;
if b > a
    bridge=fun_TanhSinh(f,a,b,tol);
    %bridge=fun_GaussKronrodBoost(f,a,b,tol);
end
%%
y=result+bridge;
end