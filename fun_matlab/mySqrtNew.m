%% SQRT Riemann Sheets
function[w]=mySqrtNew(z)
% w       =   sqrt(z);
% w       =   w.*(imag(w)<=0)-w.*(imag(w)>0); % Proper Sheet
w=sqrt(z);
for ii = 1:length(z)
if imag(w(ii))>0
    w(ii)=real(w(ii))-1j*imag(w(ii));
end
if real(w(ii))<0
   w(ii)=-real(w(ii))+1j*imag(w(ii)); 
end
end
end