%% Brief Levin transformation function with the Sidi W-algorithm extrapolation
function [val,A,B] = fun_LevinSidi(k,Sk,wk,X,A,B)
B(k-1)=1.0/wk;
A(k-1)=Sk*B(k-1);
for j=1:k-2
    d=1.0/X(k)-1.0/X(k-j);
    A(k-j-1)=(A(k-j)-A(k-j-1))/d;
    B(k-j-1)=(B(k-j)-B(k-j-1))/d;
end
val=A(1)/B(1);
end