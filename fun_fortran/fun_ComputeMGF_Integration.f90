!! MGF Integration
!Dyadic green tensor 3x3 and dyadic phi
subroutine fun_ComputeMGF_Integration(loc_r1,loc_r,loc_e,loc_k,freq,GA,Gphi)
use myfargs
implicit none
!Input
real*8 loc_r1(3), loc_r(3), r_perp(2)
real*8 freq
complex*16 loc_k(2), loc_e(2)
!Output
complex*16 GA(3,3), Gphi, near, far 
!User
real*8 swithching_point, knorm(2), knorm_max
real*8 csi, cos_term, sin_term
complex*16 KG(5)
double complex fun_IntegrateSpectralFarField, fun_IntegrateSpectralNearField
external fun_IntegrateSpectralFarField, fun_IntegrateSpectralNearField
!
r1=loc_r1
r=loc_r
k=loc_k
e=loc_e
!write(6,*) 'k(1)=', loc_k(1)
!write(6,*) 'k(2)=', loc_k(2)
KG(1:5)=cmplx(0.0,0.0) !components of Dyadic Green
! Swithching point
knorm(1)=abs(k(1))
knorm(2)=abs(k(2))
knorm_max=maxval(knorm);
swithching_point=-1.d0; !default
if (swithching_point<0) then
    swithching_point=1.2*knorm_max;
end if
! Distance perp
r_perp=r(1:2)-r1(1:2)
rho=norm2(r_perp)
!write(6,*) 'rho=', rho
!write(6,*) 'freq=', freq
! Compute Dyadic components
near = fun_IntegrateSpectralNearField(r1,r,e,k,freq,rho,1,0.d0,swithching_point)
far  =  fun_IntegrateSpectralFarField(r1,r,e,k,freq,rho,1,0.d0,swithching_point)
!write(66,*) 'K1'
!write(66,*) 'near', near
!write(66,*) 'far', far
KG(1)=near+far
!
near = fun_IntegrateSpectralNearField(r1,r,e,k,freq,rho,2,1.d0,swithching_point)
far  = fun_IntegrateSpectralFarField(r1,r,e,k,freq,rho,2,1.d0,swithching_point)
!write(66,*) 'K2'
!write(66,*) 'near', near
!write(66,*) 'far', far
KG(2)=near+far
!
near = fun_IntegrateSpectralNearField(r1,r,e,k,freq,rho,3,1.d0,swithching_point)
far  = fun_IntegrateSpectralFarField(r1,r,e,k,freq,rho,3,1.d0,swithching_point)
!write(66,*) 'K3'
!write(66,*) 'near', near
!write(66,*) 'far', far
KG(3)=near+far
!
near = fun_IntegrateSpectralNearField(r1,r,e,k,freq,rho,4,0.d0,swithching_point)
far  = fun_IntegrateSpectralFarField(r1,r,e,k,freq,rho,4,0.d0,swithching_point)
!write(66,*) 'K4'
!write(66,*) 'near', near
!write(66,*) 'far', far
KG(4)=near+far
!
far=fun_IntegrateSpectralFarField(r1,r,e,k,freq,rho,5,0.d0,swithching_point)
near=fun_IntegrateSpectralNearField(r1,r,e,k,freq,rho,5,0.d0,swithching_point)
!write(66,*) 'K5'
!write(66,*) 'near', near
!write(66,*) 'far', far
KG(5)=near+far
! Set up Tensor
csi=atan2((r(2)-r1(2)),(r(1)-r1(1)))
cos_term=cos(csi);
sin_term=sin(csi);
!
GA(1,1)=KG(1)
GA(1,2)=0
GA(1,3)=KG(2)*cos_term
GA(2,1)=0
GA(2,2)=KG(1)
GA(2,3)=KG(2)*sin_term
GA(3,1)=KG(3)*cos_term
GA(3,2)=KG(3)*sin_term
GA(3,3)=KG(4)
!
Gphi=KG(5)
end subroutine fun_ComputeMGF_Integration
