!Fourth component of dyadic tensor: GA(4)
!for points in different layers
double precision function fun_GA4_mi_Im(kr) result(y)
use myfargs
implicit none
!Input
real*8 kr
!User
complex*16 krho
complex*16 fun_Ih, fun_Ie, fun_Gmi_Ie, fun_Gmi_Ih, fun, fun2
complex*16 Tmi
complex*16 Zh1, Zh2, Ze1, Ze2, Yh1, Yh2, Ye1, Ye2
complex*16 factYe, factYh
complex*16 em, ei, km, kz1, kz2, Zei, Zem, Yei, Yem, Zhi, Zhm, kzi, kzm, ki, Yhi, Yhm
complex*16 j, argB
complex*16 segno
integer NZ, IERR
complex*16 c3, c4, c5
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
Ze1=kz1/(omega*e(1))
Ze2=kz2/(omega*e(2))
Zh1=(omega*mu0)/kz1
Zh2=(omega*mu0)/kz2
!Admittance
Ye1=(omega*e(1))/kz1
Ye2=(omega*e(2))/kz2
Yh1=kz1/(omega*mu0)
Yh2=kz2/(omega*mu0)
!Reflection coefficient
factYh=(Yh2-Yh1)/(Yh2+Yh1)
factYe=(Ye2-Ye1)/(Ye2+Ye1)
if (i==1) then !Gmi = G21
  em=e(2)
  ei=e(1)
  km=k(2)
  ki=k(1)
  Zhi=Zh1
  Zei=Ze1
  Zhm=Zh2
  Zem=Ze2
  Yhi=Yh1
  Yhm=Yh2
  Yei=Ye1
  Yem=Ye2
  kzi=kz1
  kzm=kz2
  segno=cmplx(1.0,0.0)
else !Gmi = G12
  em=e(1)
  ei=e(2)
  km=k(1)
  ki=k(2)
  Zhi=Zh2
  Zei=Ze2
  Zhm=Zh1
  Zem=Ze1
  Yhi=Yh2
  Yhm=Yh1
  Yei=Ye2
  Yem=Ye1
  kzi=kz2
  kzm=kz1
  segno=cmplx(-1.0,0.0)
end if
!Transmission coefficient
Tmi=exp(-j*kzm*(-segno*r(3)))
!Gii(0,r1)
fun_Ih=0.5d0*Yhi*( exp(-j*kzi*abs(-r1(3)))+segno*factYh*exp(-j*kzi*(segno*r1(3))))
fun_Ie=0.5d0*Yei*( exp(-j*kzi*abs(-r1(3)))+segno*factYe*exp(-j*kzi*(segno*r1(3))))
!Gmi
fun_Gmi_Ie=fun_Ie*Tmi
fun_Gmi_Ih=fun_Ih*Tmi
!Coefficients
c3=ki**2
c4=(kzm**2)/(km**2)
c5=mu0/ei
!Compute Green
fun=(c5/(j*omega))*(fun_Gmi_Ie-(c3/(krho**2)*(c4*fun_Gmi_Ie-fun_Gmi_Ih)))/mu0
!Compute besselJ(complex argument)
argB=cmplx(rho,0.d0)*krho
call ZBESJ(real(argB), imag(argB),order, 1, 1, CYR, CYI, NZ, IERR)
!
fun2=cmplx(CYR(1),CYI(1))*fun*(krho**(order+1))/(2*pi)
y=imag(fun2)
end function fun_GA4_mi_Im
