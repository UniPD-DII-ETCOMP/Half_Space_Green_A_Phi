complex*16 function fun_IntegrateSpectralNearField(loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,component,loc_order,a) result(y)
use myfargs
implicit none
!Input
real*8, intent(in) :: loc_r1(3), loc_r(3)
complex*16, intent(in) :: loc_e(2), loc_k(2)
real*8, intent(in) ::loc_rho, loc_freq, a
integer, intent(in) :: component
real*8, intent(in) :: loc_order
!User
complex*16 first_path, second_path, third_path
real*8 h, tol, zero_offset, knorm(2), knorm_max
complex*16, external :: fun_GaussKronrodBoost
!
knorm(1)=abs(loc_k(1))
knorm(2)=abs(loc_k(2))
knorm_max=maxval(knorm);
h=1e-2*knorm_max
tol=1e-6
zero_offset=0.d0
!!Integration paths
!First path
fixRe_flag=.TRUE.
fixed_part=zero_offset
first_path=fun_GaussKronrodBoost(0.d0,h,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
fixRe_flag=.FALSE.
!Second path
fixIm_flag=.TRUE.
fixed_part=h
second_path=fun_GaussKronrodBoost(zero_offset,a,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
fixIm_flag=.FALSE.
!Third path
fixRe_flag=.TRUE.
fixed_part=a
third_path=fun_GaussKronrodBoost(h,0.d0,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
fixRe_flag=.FALSE.
!
y=first_path+second_path+third_path
end function fun_IntegrateSpectralNearField
