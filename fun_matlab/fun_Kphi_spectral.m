function [y] = fun_Kphi_spectral(kr,r1,r,k,freq,e,mu0)
%% Definition
omega=2*pi*freq;
eps0=8.8541878128e-12;
%% Wave number
kz1=@(kr) mySqrtNew(k(1).^2-kr.^2); %source media
kz2=@(kr) mySqrtNew(k(2).^2-kr.^2); %target media
%% Impedenza
Ze1=@(kr) kz1(kr)./(omega*e(1));
Ze2=@(kr) kz2(kr)./(omega*e(2));

Zh1=@(kr) (omega*mu0)./kz1(kr);
Zh2=@(kr) (omega*mu0)./kz2(kr);

factZe=@(kr) (Ze2(kr)-Ze1(kr))./(Ze2(kr)+Ze1(kr));
factZh=@(kr) (Zh2(kr)-Zh1(kr))./(Zh2(kr)+Zh1(kr));
%% 
zzam=abs(r(3)-r1(3));
zzp=(r(3)+r1(3));
%%
if (r1(3)>0 && r(3)>0) %G11
    %% G11
    sGmi_Ve = @(kr) 0.5*Ze1(kr).*(...
             exp(-1j*kz1(kr)*zzam)+factZe(kr).*exp(-1j*kz1(kr)*zzp)...
             );
         
    sGmi_Vh = @(kr) 0.5*Zh1(kr).*(...
             exp(-1j*kz1(kr)*zzam)+factZh(kr).*exp(-1j*kz1(kr)*zzp)...
             );                  
elseif (r1(3)<0 && r(3)<=0) %G22
    %% G22
    sGmi_Ve = @(kr) 0.5*Ze2(kr).*(...
             exp(-1j*kz2(kr)*zzam)-factZe(kr).*exp(-1j*kz2(kr)*(-zzp))...
             );
         
    sGmi_Vh = @(kr) 0.5*Zh2(kr).*(...
             exp(-1j*kz2(kr)*zzam)-factZh(kr).*exp(-1j*kz2(kr)*(-zzp))...
             ); 
elseif (r1(3)>0 && r(3)<=0) %G21
    %% T21 coefficient
    Tmi = @(kr) exp(-1j*kz2(kr)*(-r(3)));
    %% G11(0,z')
    sGii_Ve = @(kr)0.5*Ze1(kr).*(...
             exp(-1j*kz1(kr)*abs(-r1(3)))+factZe(kr).*exp(-1j*kz1(kr)*(r1(3)))...
             );
    sGii_Vh = @(kr) 0.5*Zh1(kr).*(...
             exp(-1j*kz1(kr)*abs(-r1(3)))+factZh(kr).*exp(-1j*kz1(kr)*(r1(3)))...
             );
    %% G21
    sGmi_Ve = @(kr) sGii_Ve(kr).*Tmi(kr);
    sGmi_Vh = @(kr) sGii_Vh(kr).*Tmi(kr);
elseif (r1(3)<0 && r(3)>=0) %G12
    %% T12 coefficient
    Tmi = @(kr) exp(-1j*kz1(kr)*(r(3)));
    %% G22(0,z')
    sGii_Ve = @(kr)0.5*Ze2(kr).*(...
             exp(-1j*kz2(kr)*abs(-r1(3)))-factZe(kr).*exp(-1j*kz2(kr)*(-r1(3)))...
             );
    sGii_Vh = @(kr) 0.5*Zh2(kr).*(...
             exp(-1j*kz2(kr)*abs(-r1(3)))-factZh(kr).*exp(-1j*kz2(kr)*(-r1(3)))...
             );
    %% G12
    sGmi_Ve = @(kr) sGii_Ve(kr).*Tmi(kr);
    sGmi_Vh = @(kr) sGii_Vh(kr).*Tmi(kr);
end
%% Function evaluation
y = (1./(kr.^2)).*(sGmi_Ve(kr)-sGmi_Vh(kr));
y = 1j*omega*eps0*y;
%y = kz1(kr).*y;
end

