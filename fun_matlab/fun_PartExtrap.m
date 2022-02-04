%% Partition extraction algorithm calling Levin-Sidi rule
function [val] = fun_PartExtrap(f,a,q,tol)
%% Initialization
kmax=10; 
eps=1e-6;
old=1e32; 
u=0; w=0; s=0;
A=zeros(kmax+1,1); B=zeros(kmax+1,1);
X=zeros(kmax+2,1);
X(1)=a;
%% Begin extrapolated integration over partitions
for kk=2:kmax+2
   X(kk)=X(kk-1)+q;
   u=fun_TanhSinh(f,X(kk-1),X(kk),eps);
%    if abs(u)==0 a che serve?!
%    end
   %% Execute extrapolation
   s=s+u;
   w=u; %Type of levin transformation
   [val,A,B]=fun_LevinSidi(kk,s,w,X,A,B);
   if kk>3 && abs(val-old) < tol*abs(val)
       break;
   end
   old=val;
end
%check with quadk
% val
% sum=0;
% for ii=2:kk
%     [result(ii-1)]=fun_GaussKronrodBoost(f,X(kk-1),X(kk),eps);
%     sum=sum+result(ii-1);
% end
% sum
if isnan(val) || isinf(val)
 val=fun_GaussKronrodBoost(f,X(1),X(kk),1e-9);
end
% q = integral(f,X(1),X(kk))
%%
% rr=linspace(X(1),X(kk),300);
% y=abs(f(rr));
% figure
% semilogy(rr,abs(f(rr)),'-xr','MarkerSize',10)
% hold on
% for ii=1:kk
%     xline(X(ii))
% end
% drawnow
% 
% 
% rr=linspace(0,X(kk),3000);
% y=abs(f(rr));
% figure
% semilogy(rr,abs(f(rr)),'-xr','MarkerSize',10)
% hold on
% for ii=1:kk
%     xline(X(ii))
% end

end