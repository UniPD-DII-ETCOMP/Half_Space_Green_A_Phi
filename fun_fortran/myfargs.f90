module myfargs
 implicit none
 real*8 r(3), r1(3), zzp, zzam, rho
 real*8 pi, omega
 integer i, m 
 real*8 order
 complex*16 k(2), e(2)
 logical fixRe_flag, fixIm_flag
 real*8 fixed_part
!!$OMP threadprivate(omega,pi, r, r1,i,m,k, e,zzp, zzam, rho,order,fixRe_flag,fixIm_flag,fixed_part)
end module
