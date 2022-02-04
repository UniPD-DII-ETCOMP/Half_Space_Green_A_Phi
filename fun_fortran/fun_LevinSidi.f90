!!Brief Levin transformation function with the Sidi W-algorithm extrapolation
subroutine fun_LevinSidi(k,Sk,wk,X,kmax,A,B,y)
!Input
integer k, kmax
real*8 X(kmax+2)
complex*16 A(kmax+1), B(kmax+1)
complex*16 Sk, wk
complex*16 y
!User
integer j
real*8 d
!
B(k-1)=1.0/wk;
A(k-1)=Sk*B(k-1);
do j=1,k-2
    d=1.0/X(k)-1.0/X(k-j);
    A(k-j-1)=(A(k-j)-A(k-j-1))/d;
    B(k-j-1)=(B(k-j)-B(k-j-1))/d;
enddo
y=A(1)/B(1);
end subroutine fun_LevinSidi
