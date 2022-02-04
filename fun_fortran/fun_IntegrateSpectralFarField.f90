complex*16 function fun_IntegrateSpectralFarField(loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,component,loc_order,a) result(y)
use myfargs
implicit none
!Input
real*8, intent(in) :: loc_r1(3), loc_r(3)
complex*16, intent(in) :: loc_e(2), loc_k(2)
real*8, intent(in) :: loc_rho, loc_freq, a
integer, intent(in) :: component
real*8, intent(in) :: loc_order
!User
integer k_param
real*8 sng, tol, val, b, q, xx
complex*16 resu, bridge
complex*16, external :: fun_GaussKronrodBoost, fun_PartExtrap, fun_TanhSinh
!
k_param=10
xx=0
sng=-log(xx)
pi=4.d0*atan(1.d0)
!
tol=1e-6
val=a*loc_rho
call fun_BesselJ_NextZero(val,loc_order,b)
b=b/loc_rho
q=pi/loc_rho
! MGF Integrand function definition
! For small rho, Bessel function oscilllations are slow, so use quadrature
if (loc_rho<1e-15) then
  resu=fun_GaussKronrodBoost(a,sng,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
else
  resu=fun_PartExtrap(b,q,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order,k_param)
end if
bridge=cmplx(0.0,0.0)
if (b>a) then
  bridge=fun_TanhSinh(a,b,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
end if
y=resu+bridge
end function fun_IntegrateSpectralFarField
