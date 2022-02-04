!!Brief Helper function for tanh-sinh quadrature
complex*16 function fun_Term(a,b,component,ekh,sigma,eta,loc_r1,loc_r,loc_e,loc_k,loc_freq,loc_rho,loc_order) result(y)
use myfargs
implicit none
!Input
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
integer component
real*8 loc_order
real*8 loc_r1(3), loc_r(3)
complex*16 loc_e(2), loc_k(2)
real*8 loc_freq, loc_rho
real*8 ekh, sigma, eta, a, b
!User
real*8 q, delta, w, freq
real*8 val1Re, val2Re, val1Im, val2Im
real*8 temp1, temp2
complex*16 val1, val2 
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
q=exp(-eta*(ekh-1.0/ekh));
delta=2.0*q/(1.0+q);
w=eta*(ekh+1.0/ekh)*delta/(1.0+q);
!
temp1=a+sigma*delta
temp2=b-sigma*delta
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
        val1Re=fun_Gphi_ii_Re(temp1)
        val1Im=fun_Gphi_ii_Im(temp1)
        val2Re=fun_Gphi_ii_Re(temp2)
        val2Im=fun_Gphi_ii_Im(temp2)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        val1Re=fun_Gphi_mi_Re(temp1)
        val1Im=fun_Gphi_mi_Im(temp1)
        val2Re=fun_Gphi_mi_Re(temp2)
        val2Im=fun_Gphi_mi_Im(temp2)
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
        val1Re=fun_GA1_ii_Re(temp1)
        val1Im=fun_GA1_ii_Im(temp1)
        val2Re=fun_GA1_ii_Re(temp2)
        val2Im=fun_GA1_ii_Im(temp2)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        val1Re=fun_GA1_mi_Re(temp1)
        val1Im=fun_GA1_mi_Im(temp1)
        val2Re=fun_GA1_mi_Re(temp2)
        val2Im=fun_GA1_mi_Im(temp2)
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
        val1Re=fun_GA2_ii_Re(temp1)
        val1Im=fun_GA2_ii_Im(temp1)
        val2Re=fun_GA2_ii_Re(temp2)
        val2Im=fun_GA2_ii_Im(temp2)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        val1Re=fun_GA2_mi_Re(temp1)
        val1Im=fun_GA2_mi_Im(temp1)
        val2Re=fun_GA2_mi_Re(temp2)
        val2Im=fun_GA2_mi_Im(temp2)
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
        val1Re=fun_GA3_ii_Re(temp1)
        val1Im=fun_GA3_ii_Im(temp1)
        val2Re=fun_GA3_ii_Re(temp2)
        val2Im=fun_GA3_ii_Im(temp2)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        val1Re=fun_GA3_mi_Re(temp1)
        val1Im=fun_GA3_mi_Im(temp1)
        val2Re=fun_GA3_mi_Re(temp2)
        val2Im=fun_GA3_mi_Im(temp2)
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
        val1Re=fun_GA4_ii_Re(temp1)
        val1Im=fun_GA4_ii_Im(temp1)
        val2Re=fun_GA4_ii_Re(temp2)
        val2Im=fun_GA4_ii_Im(temp2)
    else if ((r1(3)>0.d0 .AND. r(3)<=0.d0) .OR. (r1(3)<0.d0 .AND. r(3)>=0.d0)) then
        if (r1(3)>0 .AND. r(3)<=0) then
            i=1
            m=2
        else if (r1(3)<0 .AND. r(3)>=0) then
            i=2
            m=1
        end if
        val1Re=fun_GA4_mi_Re(temp1)
        val1Im=fun_GA4_mi_Im(temp1)
        val2Re=fun_GA4_mi_Re(temp2)
        val2Im=fun_GA4_mi_Im(temp2)
    end if
end if
!
val1=cmplx(val1Re,val1Im)
val2=cmplx(val2Re,val2Im)
y=w*(val1+val2)
end function fun_Term
