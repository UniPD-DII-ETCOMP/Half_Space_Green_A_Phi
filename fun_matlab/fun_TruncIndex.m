%% Brief Helper function for tanh-sinh quadrature
function [n,s] = fun_TruncIndex(f,eh,s,sigma,eta,a,b)
nmax=24;
kappa=1e-15;
ekh=eh;
%%
for n=0:nmax
    t=fun_Term(f,ekh,sigma,eta,a,b);
    s=s+t;
    if abs(t) <= kappa*abs(s)
        break;
    end
    ekh=ekh*eh;
end
n=n-1;
end