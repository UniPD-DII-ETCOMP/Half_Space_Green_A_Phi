!!SQRT Riemann Sheets
subroutine mySqrtNew(z,w)
implicit none
complex*16 z,w
real*8 tmp1, tmp2
!
w=sqrt(z)
if  (imag(w) > 0.d0) then
    tmp1=real(w)
    tmp2=-imag(w)
    w=cmplx(tmp1,tmp2)
end if
if (real(w) < 0.d0) then
    tmp1=-real(w)
    tmp2=imag(w)
    w=cmplx(tmp1,tmp2)
end if
end subroutine mySqrtNew
