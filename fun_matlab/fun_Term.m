%% Brief Helper function for tanh-sinh quadrature
function [t] = fun_Term(f,ekh,sigma,eta,a,b)
q=exp(-eta*(ekh-1.0/ekh));
delta=2.0*q/(1.0+q);
w=eta*(ekh+1.0/ekh)*delta/(1.0+q);
t=w*(f(a,sigma*delta)+f(b,-sigma*delta));
end