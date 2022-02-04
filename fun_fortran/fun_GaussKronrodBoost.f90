!!Brief Wrapped function that uses Boost's Gauss-Kronrod quadrature routine
!to integrate a given complex-valued function along the real line.
complex*16 function fun_GaussKronrodBoost(loc_a,loc_b,tol,component,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order) result(y)
use myfargs
implicit none
!Input
real*8, intent(in) :: loc_a, loc_b, tol
integer, intent(in) :: component
real*8, intent(in) :: loc_order
real*8, intent(in) :: loc_r1(3), loc_r(3)
complex*16, intent(in) :: loc_e(2), loc_k(2)
real*8, intent(in) :: loc_freq, loc_rho
!User
real*8 freq, c, a, b
integer neval, ier, limit, lenw,last
real*8 abserr,epsabs,epsrel
real*8 error, multiplier, resuRe, resuIm
integer,allocatable,dimension(:) :: iwork
double precision,allocatable,dimension(:) :: work
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
!integer iwork(1000)
!real*8 work(5000)
!
rho=loc_rho
freq=loc_freq
r=loc_r
r1=loc_r1
k=loc_k
e=loc_e
order=loc_order
pi=4.d0*atan(1.d0)
omega=2*pi*freq
!
error=1.d0
epsabs=1e-5
epsrel=tol
limit=1000
multiplier=1.d0
if (loc_a > loc_b) then
  c=loc_b
  b=loc_a
  a=c
  multiplier=-1.d0
else 
    a=loc_a
    b=loc_b
end if
!Integrator settings
allocate(iwork(limit))
lenw=limit*5
allocate(work(lenw))
!Selection of integrand function
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
        call dqags(fun_Gphi_ii_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_Gphi_ii_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        call dqags(fun_Gphi_mi_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_Gphi_mi_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
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
        call dqags(fun_GA1_ii_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_GA1_ii_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        call dqags(fun_GA1_mi_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_GA1_mi_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
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
        call dqags(fun_GA2_ii_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_GA2_ii_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        call dqags(fun_GA2_mi_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_GA2_mi_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
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
        call dqags(fun_GA3_ii_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_GA3_ii_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        call dqags(fun_GA3_mi_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_GA3_mi_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
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
        call dqags(fun_GA4_ii_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_GA4_ii_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        call dqags(fun_GA4_mi_Re,a,b,epsabs,epsrel,resuRe,abserr,neval,ier,limit,lenw,last,iwork,work)
        call dqags(fun_GA4_mi_Im,a,b,epsabs,epsrel,resuIm,abserr,neval,ier,limit,lenw,last,iwork,work)
        y=multiplier*cmplx(resuRe,resuIm)
    end if
end if
end function fun_GaussKronrodBoost
