!!Brief Helper function for tanh-sinh quadrature
subroutine fun_TruncIndex(a,b,component,eh,loc_s,sigma,eta,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order,nn,ss)
use myfargs
implicit none
!Input
integer, intent(in) :: component
real*8, intent(in) :: loc_order
real*8, intent(in) :: loc_r1(3), loc_r(3)
complex*16, intent(in) :: loc_e(2), loc_k(2)
real*8, intent(in) :: loc_freq, loc_rho
complex*16, intent(in) :: loc_s
real*8, intent(in) :: eh, sigma, eta, a, b
integer nn
complex*16 ss
!User
integer nmax, n
real*8 kappa, ekh
complex*16 term, s
complex*16, external :: fun_Term
s=loc_s
!
nmax=24
kappa=1e-15
ekh=eh
do n=0,nmax
    term=fun_Term(a,b,component,ekh,sigma,eta,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
    s=s+term
    if (abs(term) <= kappa*abs(s)) then
        exit
    end if
    ekh=ekh*eh
enddo
n=n-1
!
nn=n
ss=s
end subroutine fun_TruncIndex
