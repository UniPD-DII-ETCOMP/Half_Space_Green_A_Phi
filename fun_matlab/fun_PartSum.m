%% Brief Helper function for tanh-sinh quadrature
function [s] = fun_PartSum(f,eh,e2h,n,sigma,eta,a,b)
ekh=eh;
s=fun_Term(f,ekh,sigma,eta,a,b);
%%
for k=2:n
   ekh=ekh*e2h;
   s=s+fun_Term(f,ekh,sigma,eta,a,b);
end
end