%% Brief Wrapped function that uses Boost's Gauss-Kronrod quadrature routine 
%to integrate a given complex-valued function along the real line. 
function [result] = fun_GaussKronrodBoost(f,a,b,tol)
%% Initialization
% order=31;
max_levels=1000;
% error=1.0;
multiplier=1.0;
%% Make sure the integration path is from the smaller number to the larger one, for Booost compatibility
if a > b
   c=b;
   b=a;
   a=c;
   multiplier=-1.0;
end
% result=integrate(f,a,b,max_levels,tol,error); %sostiruire con quadk
%,'MaxIntervalCount',max_levels

[result,errbnd]=quadgk(f,a,b,'RelTol',tol,'MaxIntervalCount',max_levels); %non sicuro 

% f_re=@(x) real(f(x));
% f_im=@(x) imag(f(x));
% 
% [result_r,errbnd]=quadgk(f_re,a,b,'RelTol',tol,'MaxIntervalCount',max_levels); %non sicuro 
% [result_i,errbnd]=quadgk(f_im,a,b,'RelTol',tol,'MaxIntervalCount',max_levels); %non sicuro 
% result=result_r+1j*result_i;
result=result*multiplier;
end