!Second component of dyadic array: GA(2)
!for points in same layer
double precision function fun_GA2_ii_Re(kr) result(y)
use myfargs
implicit none
!Input
real*8 kr
!User
complex*16 krho
complex*16 fun_Vh, fun_Ve, fun, fun2
complex*16 kz1, kz2, Ze1, Ze2, Zh1, Zh2, factZe, factZh, Zei, Zhi, kzi
complex*16 ki, ei, c1, c2
complex*16 j, argB
complex*16 segno
integer NZ, IERR
real*8 mu0
double precision CYR(1), CYI(1)
!
krho=cmplx(kr,0.0)
!
j=cmplx(0.d0,1.d0)
pi=4.d0*atan(1.d0)
mu0=4.d0*pi*1e-7
!
zzam=abs(r(3)-r1(3))
zzp=(r(3)+r1(3))
!Modify kr if required
if (fixRe_flag .OR. fixIm_flag) then
  if (fixRe_flag) then
    krho=cmplx(fixed_part,kr)
  else if (fixIm_flag) then
    krho=cmplx(kr,fixed_part)
  end if
end if
!Wave number
call mySqrtNew(k(1)**2.d0-krho**2.d0,kz1)
call mySqrtNew(k(2)**2.d0-krho**2.d0,kz2)
!Impedance
Ze1=kz1/(omega*e(1))
Ze2=kz2/(omega*e(2))
Zh1=(omega*mu0)/kz1
Zh2=(omega*mu0)/kz2
!Reflection coefficient
factZe=(Ze2-Ze1)/(Ze2+Ze1)
factZh=(Zh2-Zh1)/(Zh2+Zh1)
! G11 or G22
if (i==1) then
  Zei=Ze1
  Zhi=Zh1
  kzi=kz1
  ki=k(1)
  ei=e(1)
  segno=cmplx(1.0,0.0)
else
  Zei=Ze2
  Zhi=Zh2
  kzi=kz2
  ki=k(2)
  ei=e(2)
  segno=cmplx(-1.0,0.0)
end if
!Compute Green
fun_Ve=segno*0.5d0*j*kzi/Zei*segno*factZe*exp(-j*kzi*(segno*zzp))
fun_Vh=segno*0.5d0*j*kzi/Zhi*segno*factZh*exp(-j*kzi*(segno*zzp))
!Constant
c1=mu0/ei
c2=(ki**2)/(kzi**2)
fun=-(c1/(j*omega))*(1/(krho**2))*(fun_Ve-c2*fun_Vh)/mu0
!Compute besselJ(complex argument)
argB=cmplx(rho,0.d0)*krho
call ZBESJ(real(argB), imag(argB),order, 1, 1, CYR, CYI, NZ, IERR)
!
fun2=cmplx(CYR(1),CYI(1))*fun*(krho**(order+1))/(2*pi)
y=real(fun2)
end function fun_GA2_ii_Re
