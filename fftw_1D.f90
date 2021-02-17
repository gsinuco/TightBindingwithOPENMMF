!  Filename : fftw_1D.f90
!  Title : Fortran program to evaluate the FFT of an array.
!           it requires the FFTW library and the file fftw3.f
!
!  Author : German Sinuco
!  Date   : 18 November 2020
!  University of Durham
!  University of Sussex

subroutine fftw_1D(n,in,out,direction,info)


  ! n    (in)     : size of the arrays in and out
  ! in   (in)     : array to evaluate the FFT
  ! out  (out)    : FFT(direction) of the array in
  !                 here we shift the elements and move th
  !                 value corresonding to the zero frequency
  !                 at the centre of the array 
  ! direction (in): 'FORWARD' or 'BACKWARD'
  ! info  (inout) : error flag


  implicit none
  
  integer,                intent(in)    :: n
  integer,                intent(inout) :: info
  character(len=*),       intent(in)    :: direction 
  complex*16,dimension(n),intent(in)    :: in
  complex*16,dimension(n),intent(out)   :: out
  
  include "fftw3.f"
  
  integer i,src,dst
  complex*16 :: tmp
  integer*8 plan_backward
  integer*8 plan_forward

  complex*16,dimension(n)  :: out_tmp
  
!  write ( *, * ) '# ONE DIMENSIONAL FFTW'
!  write ( *, * ) '#', direction  

  SELECT CASE (direction)
  CASE("FORWARD")
     
     !write(*,905)(real(in(I)),aimag(in(I)),I=1,N)
     call dfftw_plan_dft_1d ( plan_forward, n, in, out, &
          &  FFTW_FORWARD, FFTW_ESTIMATE )  
     call dfftw_execute ( plan_forward )
     
     call dfftw_destroy_plan ( plan_forward )

     ! RUDIMENTARY FFTSHIFT
     ! WHICH MOVE THE VALUE OF THE ZERO FREQUENCY 
     ! TO THE CENTER OF THE OUT ARRAY
     out_tmp = out
     out(n/2+1:n) = out_tmp(1:n/2)
     do i = 1,n/2
        out(i) = out_tmp(n/2+i)
     end do

  CASE("BACKWARD")
     call dfftw_plan_dft_1d ( plan_backward, n, in, out, &
          &  FFTW_BACKWARD, FFTW_ESTIMATE )
     
     call dfftw_execute ( plan_backward )
     
     call dfftw_destroy_plan ( plan_backward )
     ! RUDIMENTARY FFTSHIFT
     ! WHICH MOVE THE VALUE OF THE ZERO FREQUENCY 
     ! TO THE CENTER OF THE OUT ARRAY

     out_tmp = out
     out(n/2+1:n) = out_tmp(1:n/2)
     do i = 1,n/2
        out(i) = out_tmp(n/2+i)
     end do

  END SELECT
905 FORMAT(2E14.6E3)

END SUBROUTINE fftw_1D

