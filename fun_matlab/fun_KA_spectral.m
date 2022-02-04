function [y] = fun_KA_spectral(kr,r1,r,k,freq,e,mu0,component) %r1 source, r target
%% Definition
omega=2*pi*freq;
%% Wave number
kz1=@(kr) mySqrtNew(k(1).^2-kr.^2); %source media
kz2=@(kr) mySqrtNew(k(2).^2-kr.^2); %target media
%% Impedenza e Ammettenza
Ze1=@(kr) kz1(kr)./(omega*e(1));
Ze2=@(kr) kz2(kr)./(omega*e(2));
Ye1=@(kr) 1./Ze1(kr);
Ye2=@(kr) 1./Ze2(kr);

Zh1=@(kr) (omega*mu0)./kz1(kr);
Zh2=@(kr) (omega*mu0)./kz2(kr);
Yh1=@(kr) 1./Zh1(kr);
Yh2=@(kr) 1./Zh2(kr);

factZe=@(kr) (Ze2(kr)-Ze1(kr))./(Ze2(kr)+Ze1(kr));
factZh=@(kr) (Zh2(kr)-Zh1(kr))./(Zh2(kr)+Zh1(kr));
factYe=@(kr) (Ye2(kr)-Ye1(kr))./(Ye2(kr)+Ye1(kr));
factYh=@(kr) (Yh2(kr)-Yh1(kr))./(Yh2(kr)+Yh1(kr));
%% 
zzam=abs(r(3)-r1(3));
zzp=(r(3)+r1(3));
%%
if (r1(3)>0 && r(3)>0) %% KA 11
    %% G11
    sGmi_Vh = @(kr) 0.5*Zh1(kr).*(...
             exp(-1j*kz1(kr)*zzam)+factZh(kr).*exp(-1j*kz1(kr)*zzp)...
             );
    sGmi_Ie = @(kr) 0.5.*(Ye1(kr)).*(...
             exp(-1j*kz1(kr)*zzam)+factYe(kr).*exp(-1j*kz1(kr)*zzp)...
             );
    sGmi_Ih = @(kr) 0.5.*(Yh1(kr)).*(...
             exp(-1j*kz1(kr)*zzam)+factYh(kr).*exp(-1j*kz1(kr)*zzp)...
             );  
    %% W11
    sWmi_Ve = @(kr) 0.5*1j*kz1(kr)./(Ze1(kr)).*factZe(kr).*exp(-1j*kz1(kr)*(zzp));
    sWmi_Vh = @(kr) 0.5*1j*kz1(kr)./(Zh1(kr)).*factZh(kr).*exp(-1j*kz1(kr)*(zzp));
    sWmi_Ie = @(kr) 0.5*1j*kz1(kr)./(Ye1(kr)).*(factYe(kr)).*exp(-1j*kz1(kr)*(zzp));
    sWmi_Ih = @(kr) 0.5*1j*kz1(kr)./(Yh1(kr)).*(factYh(kr)).*exp(-1j*kz1(kr)*(zzp));   
    %% Constants
    c1 = mu0/e(1);
    c2 = @(kr) (k(1)^2)./(kz1(kr).^2);
    c3 = k(1)^2;
    c4 = @(kr) (kz1(kr).^2)./(k(1)^2);
    c5 = c1;
elseif (r1(3)<0 && r(3)<=0) %% KA 22
    %% G22
    sGmi_Vh = @(kr) 0.5.*Zh2(kr).*(...
         exp(-1j*kz2(kr)*zzam)-factZh(kr).*exp(-1j*kz2(kr)*(-zzp))...
         ); 
    sGmi_Ie = @(kr) 0.5.*(Ye2(kr)).*(...
              exp(-1j*kz2(kr)*zzam)+(-factYe(kr)).*exp(-1j*kz2(kr)*(-zzp))...
                ); 
    sGmi_Ih = @(kr) 0.5.*(Yh2(kr)).*(...
               exp(-1j*kz2(kr)*zzam)+(-factYh(kr)).*exp(-1j*kz2(kr)*(-zzp))...
                );
    %% W22
    sWmi_Ve = @(kr) -0.5*1j*kz2(kr)./(Ze2(kr)).*(-factZe(kr)).*exp(-1j*kz2(kr)*(-(zzp)));
    sWmi_Vh = @(kr) -0.5*1j*kz2(kr)./(Zh2(kr)).*(-factZh(kr)).*exp(-1j*kz2(kr)*(-(zzp)));
    sWmi_Ie = @(kr) -0.5*1j*kz2(kr)./(Ye2(kr)).*(-factYe(kr)).*exp(-1j*kz2(kr)*(-(zzp)));
    sWmi_Ih = @(kr) -0.5*1j*kz2(kr)./(Yh2(kr)).*(-factYh(kr)).*exp(-1j*kz2(kr)*(-(zzp)));
    %% Constants 
    c1 = mu0/e(2);
    c2 = @(kr) (k(2)^2)./(kz2(kr).^2);
    c3 = k(2)^2;
    c4 = @(kr) (kz2(kr).^2)./(k(2)^2);
    c5 = c1;
elseif (r1(3)>0 && r(3)<=0) %% KA 21
    %% T21 coefficient
    Tmi = @(kr) exp(-1j*kz2(kr)*(-r(3)));
    %% G11(0,z')
    sGii_Vh = @(kr) 0.5.*Zh1(kr).*(...
             exp(-1j*kz1(kr)*abs(-r1(3)))+factZh(kr).*exp(-1j*kz1(kr)*r1(3))...
             );   
    sGii_Ie = @(kr) 0.5.*(Ye1(kr)).*(...
                exp(-1j*kz1(kr)*abs(-r1(3)))+factYe(kr).*exp(-1j*kz1(kr)*r1(3))...
                );
    sGii_Ve = @(kr) 0.5.*(Ze1(kr)).*(...
                exp(-1j*kz1(kr)*abs(-r1(3)))+factZe(kr).*exp(-1j*kz1(kr)*r1(3))...
                );   
    sGii_Ih = @(kr) 0.5.*(Yh1(kr)).*(...
             exp(-1j*kz1(kr)*abs(-r1(3)))+factYh(kr).*exp(-1j*kz1(kr)*r1(3))...
             );
    %% G21
    sGmi_Vh = @(kr) sGii_Vh(kr).*Tmi(kr);
    sGmi_Ie = @(kr) sGii_Ie(kr).*Tmi(kr);
    sGmi_Ih = @(kr) sGii_Ih(kr).*Tmi(kr);
    %% I21
    sImi_Ve = @(kr) -Ze2(kr).*sGii_Ie(kr).*Tmi(kr);
    sImi_Vh = @(kr) -Zh2(kr).*sGii_Ih(kr).*Tmi(kr);
    sImi_Ie = @(kr) -Ye2(kr).*sGii_Ve(kr).*Tmi(kr);
    sImi_Ih = @(kr) -Yh2(kr).*sGii_Vh(kr).*Tmi(kr);
    %% W21
    sWmi_Ve = @(kr) -1j*kz2(kr)./Ze2(kr).*sImi_Ve(kr);
    sWmi_Vh = @(kr) -1j*kz2(kr)./Zh2(kr).*sImi_Vh(kr);
    sWmi_Ie = @(kr) -1j*kz2(kr)./Ye2(kr).*sImi_Ie(kr);
    sWmi_Ih = @(kr) -1j*kz2(kr)./Yh2(kr).*sImi_Ih(kr);
    %% Constants
    c1 = mu0/e(2);
    c2 = @(kr) (k(2)^2)./(kz2(kr).^2);
    c3 = k(1)^2;
    c4 = @(kr) (kz2(kr).^2)./(k(2)^2);
    c5 = mu0/e(1);
elseif (r1(3)<0 && r(3)>=0) %% KA 12
    %% T12 coefficient
    Tmi = @(kr) exp(-1j*kz1(kr)*(r(3)));
    %% G22(0,z')
    sGii_Vh = @(kr) 0.5.*Zh2(kr).*(...
                exp(-1j*kz2(kr)*abs(-r1(3)))-factZh(kr).*exp(-1j*kz2(kr)*(-r1(3)))...
                );
    sGii_Ve = @(kr) 0.5.*(Ze2(kr)).*(...
                exp(-1j*kz2(kr)*abs(-r1(3)))-factZe(kr).*exp(-1j*kz2(kr)*(-r1(3)))...
                );
    sGii_Ih = @(kr) 0.5.*(Yh2(kr)).*(...
                exp(-1j*kz2(kr)*abs(-r1(3)))-factYh(kr).*exp(-1j*kz2(kr)*(-r1(3)))...
                ); 
    sGii_Ie = @(kr) 0.5.*(Ye2(kr)).*(...
                exp(-1j*kz2(kr)*abs(-r1(3)))-factYe(kr).*exp(-1j*kz2(kr)*(-r1(3)))...
                );  
    %% G12
    sGmi_Vh = @(kr) sGii_Vh(kr).*Tmi(kr);
    sGmi_Ih = @(kr) sGii_Ih(kr).*Tmi(kr);
    sGmi_Ie = @(kr) sGii_Ie(kr).*Tmi(kr);
    %% I12
    sImi_Vh = @(kr) Zh1(kr).*sGii_Ih(kr).*Tmi(kr);
    sImi_Ve = @(kr) Ze1(kr).*sGii_Ie(kr).*Tmi(kr);
    sImi_Ih = @(kr) Yh1(kr).*sGii_Vh(kr).*Tmi(kr);
    sImi_Ie = @(kr) Ye1(kr).*sGii_Ve(kr).*Tmi(kr);
    %% W12
    sWmi_Vh = @(kr) -1j*kz1(kr)./Zh1(kr).*sImi_Vh(kr);
    sWmi_Ve = @(kr) -1j*kz1(kr)./Ze1(kr).*sImi_Ve(kr);
    sWmi_Ih = @(kr) -1j*kz1(kr)./Yh1(kr).*sImi_Ih(kr);
    sWmi_Ie = @(kr) -1j*kz1(kr)./Ye1(kr).*sImi_Ie(kr);
    %% Constants
    c1 = mu0/e(1);
    c2 = @(kr) (k(1)^2)./(kz1(kr).^2);
    c3 = k(2)^2;
    c4 = @(kr) (kz1(kr).^2)./(k(1)^2);    
    c5 = mu0/e(2);
end
%% Components evaluations
%% xx, yy
if component == 1
    y = (1/(1j*omega))*sGmi_Vh(kr)/mu0;
%% xz , yz
elseif component == 2
    y = -(c1/(1j*omega))*(1./(kr.^2)).*(sWmi_Ve(kr)-c2(kr).*sWmi_Vh(kr))/mu0;
%% zx , zy
elseif component == 3
    y = -(1/(1j*omega))*(1./(kr.^2)).*(c2(kr).*sWmi_Ie(kr)-sWmi_Ih(kr))/mu0;
%% zz
elseif component == 4
    y = (c5/(1j*omega))*(sGmi_Ie(kr) - (c3./(kr.^2).*(c4(kr).*sGmi_Ie(kr)-sGmi_Ih(kr))))/mu0;
end
%y = kz1(kr).*y;
end
