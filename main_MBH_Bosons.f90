! 17th February 2021
! This programm evaluates the time evolution associated
! with the effective Floquet Hamiltonian of a driven optical
! lattice in a harmonic trap.
! I want to investigate whether effects of the trap in the 
! adiabatic deformation of the ground state.
! Without interactions or trapping potential, the "normal" ground state has momentum zero.
! When the tunneling rate is negative, the ground state has momentum pi.
!  Taking as the initial 
! condition the "normal" ground state (q=0) and ramping down the tunneling
! rate, we can't see a change of state since q=0 state is also a eigenstate of the
! Hamiltonian for any value of the tunneling rate.
!
! The function TIGHT_BINDING_HAMILTONIAN, sets the DC component of the single particle Hamiltonian 
! following Trombettoni (2001,2003) and my corrections:

! H = sum_j epsilon_j a_j a_j
!    + sum_j K_j e^phi a_j^daggert a_{j+1}+e^-phi a_{j+1}^daggert a_{j}
!    
! with
! epsilon_j = int d^3r |PHI_J(r)|^2 V_{trap}(r) ! position dependent energy shift
! K_j       = k_0 + DELTA_j                     ! position dependent tunnelling rate
! K_0       =                                   ! homoeneous tunnelling rate 
! DELTA_j   = int d^3 PHI_j PHI_{j+/-1} V_{trap}(r) ! trap effect on the tunnelling
! phi       = (m omega_mod A d/hbar) *(1 + (omega_trap/omega_mod)^2) ! 
!
! H_0 =  sum_j epsilon_j a_j a_j
!    + sum_j K_j J_0(a) (a_j^daggert a_{j+1} + a_{j+1}^daggert a_{j})
!     
! and we use the time-dependent a(t) = linear ramp. The initial state is the ground state of the lattice with a(t) = 0.
! I want so see how is it possible to make the final state close to the ground state of the final hamiltonia. Maybe using
! the fidelity.



PROGRAM TIGH_BINDING

  IMPLICIT NONE
  INTEGER INFO, N_BODIES,N_SITES,N_T,D_BARE,I_,J_
  DOUBLE PRECISION TIME,D_H,DT,NORM

  
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: H
  COMPLEX*16, DIMENSION(:),   ALLOCATABLE :: PSI
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E_BARE
  complex*16, dimension(32) :: test
  INFO     = 0

  N_BODIES = 1     
  N_SITES  = 32
  D_BARE   = D_H(N_SITES,N_BODIES,'B')  ! D_H EVALUATES THE NUMBER OF STATES  

  ALLOCATE(H(D_BARE,D_BARE))
  ALLOCATE(E_BARE(D_BARE))

  H = DCMPLX(0.0,0.0)

  N_T  = 8*1024  ! number of time steps
  DT   = 1.0/N_T ! time step
  TIME = 0.0     ! initial time
  test = 0           
  DO I_=1,N_SITES
     test(I_) = exp(-1.0*(i_-16)**2/16.0)
  END DO
  !PSI = 0
  CALL TIGHT_BINDING_HAMILTONIAN(SIZE(H,1),TIME,H,INFO)   
  CALL LAPACK_FULLEIGENVALUES(H,D_BARE,E_BARE,INFO)
  !CALL WRITE_MATRIX(real(H))
  !psi      = 0.0
  PSI      = (H(:,1))! + H(:,2) + H(:,3))/sqrt(3.0)! INITIAL CONDITION: THE WAVEFUNCTION OF THE GROUND STATE
  !DO i_=1,32!size(psi,1)
  !   write(*,*) i_
  !   !PSI(i_)  = DCMPLX(1.0,0.0)
  !END DO
  NORM     = DOT_PRODUCT(PSI, PSI)
  PSI      = PSI/SQRT(NORM)
  DO J_=1,SIZE(E_BARE,1)
     WRITE(*,*) TIME,J_,REAL(PSI(J_)),AIMAG(PSI(J_))
  END DO
  !WRITE(*,*)
  !WRITE(*,*)
  DO i_=1,16!N_T
     TIME = i_*1.0/N_T
     !CALL TIGHT_BINDING_HAMILTONIAN(SIZE(H,1),TIME,H,INFO)   
     !CALL LAPACK_FULLEIGENVALUES(H,D_BARE,E_BARE,INFO)

     !DO J_=1,SIZE(E_BARE,1)
     !   WRITE(*,*) J_,E_BARE(J_)
     !END DO
     !WRITE(*,*)
     !WRITE(*,*)
     !The eigenvectors are the columns of H_BARE, i.e.
     !<n|e_i> = H_BARE(n,i)
     !where |n> is the n-th lattice site. This is because:
     ! |e_i> = SUM_n H_BARE(n,i) |n>
     !CALL WRITE_MATRIX(real(H))
     CALL RKFOURTHORDER(TIGHT_BINDING_HAMILTONIAN,SIZE(E_BARE,1),PSI,TIME,DT,INFO)
     !CALL RKFOURTHORDER(SIZE(E_BARE,1),PSI,TIME,DT,INFO)
     IF(MOD(I_,64).EQ.0) THEN
        DO J_=1,SIZE(E_BARE,1)
           WRITE(*,*) time,J_,REAL(PSI(J_)),AIMAG(PSI(J_))
        END DO
        WRITE(*,*)
        WRITE(*,*)
     END IF
 
  END DO
 
END PROGRAM TIGH_BINDING
  

