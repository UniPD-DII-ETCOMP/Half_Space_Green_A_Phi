!!Brief Helper function for tanh-sinh quadrature
complex*16 function fun_PartSum(a,b,component,eh,e2h,n,sigma,eta,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order) result(y)
use myfargs
implicit none
!Input
integer, intent(in) :: component
real*8, intent(in) ::  loc_order
real*8, intent(in) :: loc_r1(3), loc_r(3)
complex*16, intent(in) :: loc_e(2), loc_k(2)
real*8, intent(in) :: loc_freq, loc_rho
real*8, intent(in) :: eh, e2h, sigma, eta, a, b
integer, intent(in) :: n
!User
complex*16, external :: fun_Term
complex*16 res
real*8 ekh
integer kk
!
ekh=eh
y=fun_Term(a,b,component,ekh,sigma,eta,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
!
do kk=2,n
   ekh=ekh*e2h;
   res=fun_Term(a,b,component,ekh,sigma,eta,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
   y=y+res
enddo
end function fun_PartSum
