!First component of dyadic tensor: GA(1,1)
!for points in different layers
double precision function fun_GA1_mi_Im(kr) result(y)
use myfargs
implicit none
!Input
real*8 kr
!User
complex*16 krho
complex*16 fun_Vh, fun_Gmi_Vh, fun, fun2
complex*16 Tmi
complex*16 kz1, kz2, Zh1, Zh2, factZh, Zhi, kzi, kzm
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
Zh1=(omega*mu0)/kz1
Zh2=(omega*mu0)/kz2
!Reflection coefficient
factZh=(Zh2-Zh1)/(Zh2+Zh1)
if (i==1) then !Gmi = G21
  Zhi=Zh1
  kzi=kz1
  kzm=kz2
  segno=cmplx(1.0,0.0)
else !Gmi = G12
  Zhi=Zh2
  kzi=kz2
  kzm=kz1
  segno=cmplx(-1.0,0.0)
end if
!Transmission coefficient
Tmi=exp(-j*kzm*(-segno*r(3)))
!Gii(0,r1)
fun_Vh=0.5d0*Zhi*( exp(-j*kzi*abs(-r1(3)))+segno*factZh*exp(-j*kzi*(segno*r1(3))))
!Gmi
fun_Gmi_Vh=fun_Vh*Tmi
!Compute Green
fun=1/(j*omega)*fun_Gmi_Vh/mu0
!Compute besselJ(complex argument)
argB=cmplx(rho,0.d0)*krho
call ZBESJ(real(argB), imag(argB),order, 1, 1, CYR, CYI, NZ, IERR)
!
fun2=cmplx(CYR(1),CYI(1))*fun*(krho**(order+1))/(2*pi)
y=imag(fun2)
end function fun_GA1_mi_Im
