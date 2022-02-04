%% Fixed imaginary part
function [f] = fun_IntegrandFixKrhoParts(x,fixed_part,which_part,rs,rt,e,k,freq,rho,component,order);
%% Modify krho to fix real or imaginary part of krho

if strcmp(which_part,'re')
    krho=fixed_part+1j*x;
elseif strcmp(which_part,'im')
    krho=x+1j*fixed_part;
end


f=fun_SommerfeldIntegrand(rs,rt,e,k,freq,rho,krho,component,order);
end
