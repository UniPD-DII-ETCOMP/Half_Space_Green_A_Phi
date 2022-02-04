double precision function fun_Gphi_mi_Re(kr) result(y)
use myfargs
implicit none
!Input
real*8 kr
!User
complex*16 krho
complex*16 Tmi, fun_Ve, fun_Vh, fun, fun2
complex*16 kz1, kz2, Ze1, Ze2, Zh1, Zh2, factZe, factZh, Zei, Zhi, kzi, kzm
complex*16 j, argB
complex*16 segno
integer NZ, IERR
real*8 eps0, mu0
double precision CYR(1), CYI(1)
!
krho=cmplx(kr,0.0)
!
pi=4.d0*atan(1.d0)
eps0=8.8541878128e-12;
mu0=4.d0*pi*1e-7
!Modify kr if required
if (fixRe_flag .OR. fixIm_flag) then
  if (fixRe_flag) then
    krho=cmplx(fixed_part,kr)
  else if (fixIm_flag) then
    krho=cmplx(kr,fixed_part)
  end if
end if
!
j=cmplx(0.d0,1.d0)
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
! G21 or G12
if (i==1) then
  Zei=Ze1
  Zhi=Zh1
  kzi=kz1
  kzm=kz2
  segno=cmplx(1.0,0.0)
else
  Zei=Ze2
  Zhi=Zh2
  kzi=kz2
  kzm=kz1
  segno=cmplx(-1.0,0.0)
end if
!Transmission coefficient
Tmi=exp(-j*kzm*(-segno*r(3)))
!Compute Green
fun_Ve=0.5d0*Zei*(exp(-j*kzi*abs(-r1(3)))+segno*factZe*exp(-j*kzi*(segno*r1(3))))
fun_Vh=0.5d0*Zhi*(exp(-j*kzi*abs(-r1(3)))+segno*factZh*exp(-j*kzi*(segno*r1(3))))
fun_Ve=fun_Ve*Tmi
fun_Vh=fun_Vh*Tmi
fun=(fun_Ve-fun_Vh)/(krho**2)
!Compute besselJ(complex argument)
argB=cmplx(rho,0.d0)*krho
call ZBESJ(real(argB), imag(argB),order, 1, 1, CYR, CYI, NZ, IERR)
!
fun2=j*omega*eps0*cmplx(CYR(1),CYI(1))*fun*(krho**(order+1))/(2*pi)
y=real(fun2)
end function fun_Gphi_mi_Re
