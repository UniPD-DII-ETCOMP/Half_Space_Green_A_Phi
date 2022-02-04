!! Partition extraction algorithm calling Levin-Sidi rule
complex*16 function fun_PartExtrap(a,q,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order,loc_kmax) result(y)
use myfargs
implicit none
!Input
real*8, intent(in) :: a, q, tol
integer, intent(in) ::  component, loc_kmax
real*8, intent(in) ::  loc_order
real*8, intent(in) :: loc_r1(3), loc_r(3)
complex*16, intent(in) :: loc_e(2), loc_k(2)
real*8, intent(in) :: loc_freq, loc_rho
!User
complex*16 AA(loc_kmax+1), BB(loc_kmax+1)
real*8 XX(loc_kmax+2)
integer kk
real*8 eps
complex*16 old, u ,w, s
complex*16, external :: fun_TanhSinh
complex*16, external :: fun_GaussKronrodBoost
! Initialize
eps=1e-6
AA(1:loc_kmax+1)=cmplx(0.0,0.0)
BB(1:loc_kmax+1)=cmplx(0.0,0.0)
XX(1:loc_kmax+2)=0.d0
old=cmplx(1e32,0.0)
u=cmplx(0.0,0.0)
w=cmplx(0.0,0.0)
s=cmplx(0.0,0.0)
XX(1)=a
!
!Begin extrapolated integration over partitions
do kk=2,loc_kmax+2
   XX(kk)=XX(kk-1)+q
   u=fun_TanhSinh(XX(kk-1),XX(kk),eps,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
    !write(6,*) 'XK(kk)', XX(kk)
   !Execute extrapolation
   s=s+u
   w=u !Type of levin transformation
   
   call fun_LevinSidi(kk,s,w,XX,loc_kmax,AA,BB,y)
    !write(6,*) 'y', y
    !write(6,*) 'old', old
    
   
   if (kk>3 .AND. abs(y-old) < tol*abs(y)) then
       exit
   end if
   old=y
enddo


if (isnan(abs(y)) .OR. (abs(y)>1.0d12)) then
 y=fun_GaussKronrodBoost(a,XX(kk-1),1d-9,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
endif

end function fun_PartExtrap