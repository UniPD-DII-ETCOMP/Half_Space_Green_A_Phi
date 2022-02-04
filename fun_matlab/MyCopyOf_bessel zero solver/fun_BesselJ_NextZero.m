%% Bessel Next Zero
%first zeros of order Bessel after point x. Bessel function of first kind
function [y] = fun_BesselJ_NextZero(order,x)
%% Setup
kind=1; %kind of Bessel fun
k=20; %number of zeros to calculate
%% Compute zeros
x0=besselzero(order,k,kind); 
while max(x0) < x %if too few zeros are calculated
    %increase number of computed zeros
    k=k+10;
    x0=besselzero(order,k,kind);
end
%% Find first zero of Bessel fun after given x
y=x0(find(x0>x,1,'first'));
end