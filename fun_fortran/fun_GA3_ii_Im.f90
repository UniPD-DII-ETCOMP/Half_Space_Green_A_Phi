!Third component of dyadic array: GA(3)
!for points in same layer
double precision function fun_GA3_ii_Im(kr) result(y)
use myfargs
implicit none
!Input
real*8 kr
!User
complex*16 krho
complex*16 fun_Ih, fun_Ie, fun , fun2
complex*16 kz1, kz2, Ye1, Ye2, Yh1, Yh2, factYe, factYh, Yei, Yhi, kzi, ki
complex*16 j, argB
complex*16 segno
integer NZ, IERR
real*8 mu0
complex*16 c2
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
!Admittance
Ye1=(omega*e(1))/kz1
Ye2=(omega*e(2))/kz2
Yh1=kz1/(omega*mu0)
Yh2=kz2/(omega*mu0)
!Reflection coefficient
factYe=(Ye2-Ye1)/(Ye2+Ye1)
factYh=(Yh2-Yh1)/(Yh2+Yh1)
! G11 or G22
if (i==1) then
  Yei=Ye1
  Yhi=Yh1
  kzi=kz1
  ki=k(1)
  segno=cmplx(1.0,0.0)
else
  Yei=Ye2
  Yhi=Yh2
  kzi=kz2
  ki=k(2)
  segno=cmplx(-1.0,0.0)
end if
!Compute Green
fun_Ie=segno*0.5d0*j*kzi/Yei*segno*factYe*exp(-j*kzi*(segno*zzp))
fun_Ih=segno*0.5d0*j*kzi/Yhi*segno*factYh*exp(-j*kzi*(segno*zzp))
!Constant
c2=(ki**2)/(kzi**2)
fun=-(1/(j*omega))*(1./(krho**2))*(c2*fun_Ie-fun_Ih)/mu0
!Compute besselJ(complex argument)
argB=cmplx(rho,0.d0)*krho
call ZBESJ(real(argB), imag(argB),order, 1, 1, CYR, CYI, NZ, IERR)
!
fun2=cmplx(CYR(1),CYI(1))*fun*(krho**(order+1))/(2*pi)
y=imag(fun2)
end function fun_GA3_ii_Im
