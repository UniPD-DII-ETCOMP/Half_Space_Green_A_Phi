%% Brief Execute double-exponential tanh-sinh quadrature based on Michalsky, Mosig
function [val] = fun_TanhSinh(f0,a,b,eps)
%% The algorithm requires to accept two arguments to handle end-point singularites
f=@(c,d) f0(c+d);
%%
eta=1.0;
maxlev=5;
%%
sigma=(b-a)/2.0;
gamma=(b+a)/2.0;
s=eta*f(gamma,0.0);
%%
h=1.5;
eh=exp(h);
%%
[n,s]=fun_TruncIndex(f,eh,s,sigma,eta,a,b);
old=sigma*h*s;
%%
for m=1:maxlev
   e2h=eh;
   h=h/2.0;
   eh=exp(h);
   s=fun_PartSum(f,eh,e2h,n,sigma,eta,a,b);
   val=old/2.0+sigma*h*s;
   if abs(val-old) <= eps*abs(val)
       break;
   end
   old=val; 
   n=n*2;
end
end