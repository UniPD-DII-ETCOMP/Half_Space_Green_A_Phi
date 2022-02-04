%% Function Dyadic Green KA, Kphi
function [KA,Kphi] = fun_DyadicGreen(r_src,r_trg,e,k,freq)
%% Homogeneous material
% if e(1) == e(2)
%     mu0=4*pi*1e-7;
%     k0=k(1);
%     [KA]=1e-7*(exp(-1j*k0*norm(r_src-r_trg))/norm(r_src-r_trg))/mu0*eye(3)
%     [Kphi]=exp(-1j*k0*norm(r_src-r_trg))/(4*pi*norm(r_src-r_trg))
% %     return
% end
% if strcmp(integration_method,'integrate')
    [K]=fun_ComputeMGF_Integration(r_src,r_trg,e,k,freq);
% elseif strcmp(integration_method,'dcim')
%     [K]=fun_ComputeMGF_DCIM(r_src,r_trg,e,k,freq);    
% elseif strcmp(integration_method,'quasistatic')
%     [K]=fun_ComputeMGF_Quasistatic(r_src,r_trg,e,k,freq);
% end
%% cos_term and sin_term
csi=atan2((r_trg(2)-r_src(2)),(r_trg(1)-r_src(1)));
cos_term=cos(csi);
sin_term=sin(csi);
%%
% if 1 
KA_spatial(1) = K(1);
KA_spatial(2) = 0;
KA_spatial(3) = K(2)*cos_term;
KA_spatial(4) = 0;
KA_spatial(5) = K(1);
KA_spatial(6) = K(2)*sin_term;
KA_spatial(7) = K(3)*cos_term;
KA_spatial(8) = K(3)*sin_term;
KA_spatial(9) = K(4);
% end
% if 0
% KA_spatial(1) = K(1);
% KA_spatial(2) = 0;
% KA_spatial(3) = K(2);
% KA_spatial(4) = 0;
% KA_spatial(5) = K(1);
% KA_spatial(6) = K(2);
% KA_spatial(7) = K(3);
% KA_spatial(8) = K(3);
% KA_spatial(9) = K(4);    
% end
%%
% KA = reshape(KA_spatial,3,3);
KA(1,1)=KA_spatial(1);
KA(1,2)=KA_spatial(2);
KA(1,3)=KA_spatial(3);
KA(2,1)=KA_spatial(4);
KA(2,2)=KA_spatial(5);
KA(2,3)=KA_spatial(6);
KA(3,1)=KA_spatial(7);
KA(3,2)=KA_spatial(8);
KA(3,3)=KA_spatial(9);
%
Kphi = K(5);
end