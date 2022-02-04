%% Test
%find first zero of first kind Bessel function after given x
clc
clear
close all
%%
x=10; %point of interest
nu=1; %order of Bessel
kind=1; %kind of Bessel fun
k=100; %number of zeros to calculate
%% Generate Bessel
z = 0:0.1:30;
J=besselj(nu,z);
%% Compute zeros
x0=besselzero(nu,k,kind);
return
%% Find zero 
% xx=findzero(nu,k,x,kind)
%% Find first zero of Bessel fun after given x
x0_first_next=x0(find(x0>x,1,'first'));
%% plot
plot(z,J)
grid on
hold on
scatter(x0,zeros(k,1),'filled','r')
scatter(x,0,'filled','b')
legend(['J_',num2str(nu)],'zeros','point of interest')
xlabel('z','interpreter','latex')
ylabel('$J_\nu(z)$','interpreter','latex')