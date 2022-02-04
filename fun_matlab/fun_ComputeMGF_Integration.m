%% MGF Integratior
function [K]=fun_ComputeMGF_Integration(rs,rt,e,k,freq)
%% Swithching point
k_max=max(abs(k));
swithching_point=-1; %default
if swithching_point<0
    swithching_point=1.2*abs(k_max);
end
%% Distance perp
rho=norm(rt(1:2)-rs(1:2));
%%
K=zeros(5,1); %all components of Dyadic Green
%%
% disp('------------')
% disp('K1')
near = fun_IntegrateSpectralNearField(rs,rt,e,k,freq,rho,1,0,swithching_point);
far  = fun_IntegrateSpectralFarField(rs,rt,e,k,freq,rho,1,0,swithching_point);
K(1)=near+far;
%% 
% disp('------------')
% disp('K2')
near = fun_IntegrateSpectralNearField(rs,rt,e,k,freq,rho,2,1,swithching_point);
far  = fun_IntegrateSpectralFarField(rs,rt,e,k,freq,rho,2,1,swithching_point);
K(2)=near+far;
%% 
% disp('------------')
% disp('K3')
near = fun_IntegrateSpectralNearField(rs,rt,e,k,freq,rho,3,1,swithching_point);
far  = fun_IntegrateSpectralFarField(rs,rt,e,k,freq,rho,3,1,swithching_point);
K(3)=near+far;
%% 
% disp('------------')
% disp('K4')
near = fun_IntegrateSpectralNearField(rs,rt,e,k,freq,rho,4,0,swithching_point);
far  = fun_IntegrateSpectralFarField(rs,rt,e,k,freq,rho,4,0,swithching_point);
K(4)=near+far;
%%
% disp('------------')
% disp('K5')
near = fun_IntegrateSpectralNearField(rs,rt,e,k,freq,rho,5,0,swithching_point);
far  = fun_IntegrateSpectralFarField(rs,rt,e,k,freq,rho,5,0,swithching_point);
K(5)=near+far;
end