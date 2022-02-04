!---- never change ------

#include "fintrf.h"      

subroutine mexFunction(nlhs, plhs, nrhs, prhs)

    implicit none

! mwPointer mexFunction arguments:
      mwPointer plhs(*), prhs(*)
      integer nlhs, nrhs
!---- end never change ------
	  
	  
! mwSize stuff for mexing
      mwSize mo,no,siz
	  
! mwPointer stuff for mexing
      mwPointer mxCreateDoubleMatrix
      mwPointer mxGetPr
      mwPointer r1_pr,r_pr,eRe_pr,eIm_pr
	  mwPointer kRe_pr,kIm_pr,freq_pr
	  mwPointer GARe_pr, GphiRe_pr
	  mwPointer GAIm_pr, GphiIm_pr
      mwPointer m, n
      mwPointer mxGetM, mxGetN

!integer (normal not 4/8)
      integer flag,i,j,ii,jj,icount
      integer mxIsNumeric 
      integer*4 ComplexFlag
	  
! fortran subroutine arguments
	  !real*8,allocatable,dimension(:,:) :: NN, G_r, Lre, Lim
      !integer,allocatable,dimension(:,:) :: G
	  !real*8,allocatable,dimension(:) :: radius 
	  real*8 r1(3),r(3),freq
	  real*8 eRe(2), eIm(2)
	  real*8 kRe(2), kIm(2)
	  complex*16 k(2), e(2)
	  real*8 GARe(3,3), GAIm(3,3)
	  real*8 GphiRe, GphiIm
	  complex*16 GA(3,3), Gphi
	  character*80 msg
      logical debu
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!      
      debu = .true. ! .true. o .false. per attivare o disattivare il debug
	  if(debu) open(unit=66,file='log.txt',status='unknown')
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    Check to see INPUTS are numeric.
	  do ii = 1,7
        if (mxIsNumeric(prhs(ii)) .ne. 1) then
          call mexErrMsgIdAndTxt ('MATLAB:fun_compute_Matrix_R:NonNumeric', &
                                'Inputs must be numeric.')
        endif
	  enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #1 is REAL and fetch it
      m = mxGetM(prhs(1))
      n = mxGetN(prhs(1))
      if(m .ne. 1 .or. n .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 1 must be 1x3.')
      endif	  
      siz = m*n
      r1_pr = mxGetPr(prhs(1))
      call mxCopyPtrToReal8(r1_pr, r1, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 1'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #2 is REAL and fetch it
      m = mxGetM(prhs(2))
      n = mxGetN(prhs(2))
      if(m .ne. 1 .or. n .ne. 3) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 2 must be scalar.')
      endif	  
      siz = m*n
      r_pr = mxGetPr(prhs(2))
      call mxCopyPtrToReal8(r_pr, r, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 2'	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #3 is REAL and fetch it
      m = mxGetM(prhs(3))
      n = mxGetN(prhs(3))
      if(m .ne. 1 .or. n .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 3 must be 1x2.')
      endif	  
      siz = m*n
      eRe_pr = mxGetPr(prhs(3))
      call mxCopyPtrToReal8(eRe_pr, eRe, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 3'	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #4 is REAL and fetch it
      m = mxGetM(prhs(4))
      n = mxGetN(prhs(4))
      if(m .ne. 1 .or. n .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 4 must be 1x2.')
      endif	  
      siz = m*n
      eIm_pr = mxGetPr(prhs(4))
      call mxCopyPtrToReal8(eIm_pr, eIm, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 4'	  	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #5 is REAL and fetch it
      m = mxGetM(prhs(5))
      n = mxGetN(prhs(5))
      if(m .ne. 1 .or. n .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 5 must be 1x2.')
      endif	  
      siz = m*n
      kRe_pr = mxGetPr(prhs(5))
      call mxCopyPtrToReal8(kRe_pr, kRe, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 5'	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	  
!     Check that input #6 is REAL and fetch it
      m = mxGetM(prhs(6))
      n = mxGetN(prhs(6))
      if(m .ne. 1 .or. n .ne. 2) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 6 must be 1x2.')
      endif	  
      siz = m*n
      kIm_pr = mxGetPr(prhs(6))
      call mxCopyPtrToReal8(kIm_pr, kIm, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 6'	  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	
 if(debu) write(66,*) 'eRe', eRe  
 if(debu) write(66,*) 'eIm', eIm  
 if(debu) write(66,*) 'kRe', kRe  
 if(debu) write(66,*) 'kIm', kIm  
	e(1)=cmplx(eRe(1),eIm(1))
	 if(debu) write(66,*) 'here0 - 1'  
	e(2)=cmplx(eRe(2),eIm(2))
	 if(debu) write(66,*) 'here0 - 2'  
    k(1)=cmplx(kRe(1),kIm(1))
	 if(debu) write(66,*) 'here0 - 3'  
	k(2)=cmplx(kRe(2),kIm(2))
     if(debu) write(66,*) 'here0 - 4'  	
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!		  
!     Check that input #7 is REAL and fetch it
      m = mxGetM(prhs(7))
      n = mxGetN(prhs(7))
      if(m .ne. 1 .or. n .ne. 1) then
         call mexErrMsgIdAndTxt ('MATLAB:convec:NonRowVector', &
                                'Input 7 must be 1x1.')
      endif	  
      siz = m*n
      freq_pr = mxGetPr(prhs(7))
      call mxCopyPtrToReal8(freq_pr, freq, siz) ! da double precision a reale
	  if(debu) write(66,*) 'input 7'		  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
! call the computational subroutine.
      call fun_ComputeMGF_Integration(r1,r,e,k,freq,GA,Gphi)

GARe=real(GA)
GAIm=imag(GA)
GphiRe=real(Gphi)
GphiIm=imag(Gphi)

       if(debu) write(66,*) 'here'
	   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	




! Create a matrix for the return argument 1
      mo=3
      no=3
	  ComplexFlag = 0
      plhs(1) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output 1 into a MATLAB array.
      GARe_pr = mxGetPr(plhs(1))
      siz=mo*no
      call mxCopyReal8ToPtr(GARe, GARe_pr, siz)
	  
	   if(debu) write(66,*) 'here2 - 1'
	  
! Create a matrix for the return argument 2
      mo=3
      no=3
	  ComplexFlag = 0
      plhs(2) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output 2 into a MATLAB array.
      GAIm_pr = mxGetPr(plhs(2))
      siz=mo*no
      call mxCopyReal8ToPtr(GAIm, GAIm_pr, siz)	  
	  
	  
	   if(debu) write(66,*) 'here2 - 2'
	  
! Create a matrix for the return argument 3
      mo=1
      no=1
	  ComplexFlag = 0
      plhs(3) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output 3 into a MATLAB array.
      GphiRe_pr = mxGetPr(plhs(3))
      siz=mo*no
      call mxCopyReal8ToPtr(GphiRe, GphiRe_pr, siz)
	  
	   if(debu) write(66,*) 'here2 - 3'

! Create a matrix for the return argument 4
      mo=1
      no=1
	  ComplexFlag = 0
      plhs(4) = mxCreateDoubleMatrix(mo, no, ComplexFlag)
! Load the output 4 into a MATLAB array.
      GphiIm_pr = mxGetPr(plhs(4))
      siz=mo*no
      call mxCopyReal8ToPtr(GphiIm, GphiIm_pr, siz)	  

      if(debu) write(66,*) 'here2 - 4'
      
      if(debu) close(66)
      return
      end

