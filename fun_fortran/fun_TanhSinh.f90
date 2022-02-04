complex*16 function fun_TanhSinh(a,b,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order) result(y)
use myfargs
implicit none
!Input
real*8, intent(in) :: a, b, tol
integer, intent(in) :: component
real*8, intent(in) :: loc_order
real*8, intent(in) :: loc_r1(3), loc_r(3)
complex*16, intent(in) :: loc_e(2), loc_k(2)
real*8, intent(in) :: loc_freq, loc_rho
!User variables
real*8 temp
complex*16 s, old, vC
real*8 eta, sigma, gamma, h, eh, e2h, valRe, valIm, freq
integer mm, maxlev, n, nn
complex*16 ss
!Dyadic components
double precision fun_Gphi_ii_Re, fun_Gphi_ii_Im, fun_Gphi_mi_Re, fun_Gphi_mi_Im
double precision fun_GA1_ii_Re, fun_GA1_ii_Im, fun_GA1_mi_Re, fun_GA1_mi_Im
double precision fun_GA2_ii_Re, fun_GA2_ii_Im, fun_GA2_mi_Re, fun_GA2_mi_Im
double precision fun_GA3_ii_Re, fun_GA3_ii_Im, fun_GA3_mi_Re, fun_GA3_mi_Im
double precision fun_GA4_ii_Re, fun_GA4_ii_Im, fun_GA4_mi_Re, fun_GA4_mi_Im
external fun_Gphi_ii_Re, fun_Gphi_ii_Im, fun_Gphi_mi_Re, fun_Gphi_mi_Im
external fun_GA1_ii_Re, fun_GA1_ii_Im, fun_GA1_mi_Re, fun_GA1_mi_Im
external fun_GA2_ii_Re, fun_GA2_ii_Im, fun_GA2_mi_Re, fun_GA2_mi_Im
external fun_GA3_ii_Re, fun_GA3_ii_Im, fun_GA3_mi_Re, fun_GA3_mi_Im
external fun_GA4_ii_Re, fun_GA4_ii_Im, fun_GA4_mi_Re, fun_GA4_mi_Im
complex*16, external :: fun_PartSum
!
rho=loc_rho
freq=loc_freq
r=loc_r
r1=loc_r1
k=loc_k
e=loc_e
order=loc_order
!
pi=4.d0*atan(1.d0)
omega=2*pi*freq
!
eta=1.0d0
maxlev=5
sigma=0.5d0*(b-a)
gamma=0.5d0*(b+a)
temp=gamma+0.d0
!Evaluate green
if (component==5) then
  !check position of points in layers
      if ((r1(3)>0.d0 .AND. r(3)>0.d0) .OR. (r1(3)<0.d0 .AND. r(3)<0.d0)) then
          if (r1(3)>0.d0 .AND. r(3)>0.d0)   then
              i=1
              m=1
          else if (r1(3)<0.d0 .AND. r(3)<0.d0) then
              i=2
              m=2
          end if
          valRe=fun_Gphi_ii_Re(temp)
          valIm=fun_Gphi_ii_Im(temp)
      else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
          if (r1(3)>0 .AND. r(3)<=0) then
              i=1
              m=2
          else if (r1(3)<0 .AND. r(3)>=0) then
              i=2
              m=1
          end if
          valRe=fun_Gphi_mi_Re(temp)
          valIm=fun_Gphi_mi_Im(temp)
      end if
  else if (component==1) then
  !check position of points in layers
      if ((r1(3)>0.d0 .AND. r(3)>0.d0) .OR. (r1(3)<0.d0 .AND. r(3)<0.d0)) then
          if (r1(3)>0.d0 .AND. r(3)>0.d0)   then
              i=1
              m=1
          else if (r1(3)<0.d0 .AND. r(3)<0.d0) then
              i=2
              m=2
          end if
          valRe=fun_GA1_ii_Re(temp)
          valIm=fun_GA1_ii_Im(temp)
      else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
          if (r1(3)>0 .AND. r(3)<=0) then
              i=1
              m=2
          else if (r1(3)<0 .AND. r(3)>=0) then
              i=2
              m=1
          end if
          valRe=fun_GA1_mi_Re(temp)
          valIm=fun_GA1_mi_Im(temp)
      end if
  else if (component==2) then
  !check position of points in layers
      if ((r1(3)>0.d0 .AND. r(3)>0.d0) .OR. (r1(3)<0.d0 .AND. r(3)<0.d0)) then
          if (r1(3)>0.d0 .AND. r(3)>0.d0)   then
              i=1
              m=1
          else if (r1(3)<0.d0 .AND. r(3)<0.d0) then
              i=2
              m=2
          end if
          valRe=fun_GA2_ii_Re(temp)
          valIm=fun_GA2_ii_Im(temp)
      else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
          if (r1(3)>0 .AND. r(3)<=0) then
              i=1
              m=2
          else if (r1(3)<0 .AND. r(3)>=0) then
              i=2
              m=1
          end if
          valRe=fun_GA2_mi_Re(temp)
          valIm=fun_GA2_mi_Im(temp)
      end if
  else if (component==3) then
  !check position of points in layers
      if ((r1(3)>0.d0 .AND. r(3)>0.d0) .OR. (r1(3)<0.d0 .AND. r(3)<0.d0)) then
          if (r1(3)>0.d0 .AND. r(3)>0.d0)   then
              i=1
              m=1
          else if (r1(3)<0.d0 .AND. r(3)<0.d0) then
              i=2
              m=2
          end if
          valRe=fun_GA3_ii_Re(temp)
          valIm=fun_GA3_ii_Im(temp)
      else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
          if (r1(3)>0 .AND. r(3)<=0) then
              i=1
              m=2
          else if (r1(3)<0 .AND. r(3)>=0) then
              i=2
              m=1
          end if
          valRe=fun_GA3_mi_Re(temp)
          valIm=fun_GA3_mi_Im(temp)
      end if
  else if (component==4) then
  !check position of points in layers
      if ((r1(3)>0.d0 .AND. r(3)>0.d0) .OR. (r1(3)<0.d0 .AND. r(3)<0.d0)) then
          if (r1(3)>0.d0 .AND. r(3)>0.d0)   then
              i=1
              m=1
          else if (r1(3)<0.d0 .AND. r(3)<0.d0) then
              i=2
              m=2
          end if
          valRe=fun_GA4_ii_Re(temp)
          valIm=fun_GA4_ii_Im(temp)
      else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
          if (r1(3)>0 .AND. r(3)<=0) then
              i=1
              m=2
          else if (r1(3)<0 .AND. r(3)>=0) then
              i=2
              m=1
          end if
          valRe=fun_GA4_mi_Re(temp)
          valIm=fun_GA4_mi_Im(temp)
      end if
  end if
vC=cmplx(valRe,valIm)
s=eta*vC
h=0.5d0
eh=exp(h)
!
call fun_TruncIndex(a,b,component,eh,s,sigma,eta,r1,r,e,k,freq,rho,loc_order,nn,ss)
n=nn
s=ss
old=sigma*h*s
!
do mm=1,maxlev
  e2h=eh
  h=h*0.5d0
  eh=exp(h)
  s=fun_PartSum(a,b,component,eh,e2h,n,sigma,eta,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order)
  y=old*0.5d0+sigma*h*s
  if (abs(y-old) <= tol*abs(y)) then
       exit
  end if
  old=y
  n=n*2
enddo
end function fun_TanhSinh
