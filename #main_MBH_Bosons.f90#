! 17th February 2021
! This programm evaluates the time evolution associated
! with the effective Floquet Hamiltonian of a driven optical
! lattice in a harmonic trap.
! I want to investigate the effects of the trap in the 
! adiabatic deformation of the ground state.
! With or witout interactions or trapping potential, the "normal" ground state has momentum zero.
! (see Creefield PRA 2008)
! When the tunneling rate is negative, the ground state has momentum pi.
! Nglecting inteeractions and tTaking as the initial 
! condition the "normal" ground state (momentum=0) and ramping down the tunneling
! rate, we can't see a change of state since q=0 state is also a eigenstate of the
! Hamiltonian for any value of the tunneling rate.
!
! The function TIGHT_BINDING_HAMILTONIAN, sets the static component of the single particle Hamiltonian 
! following Trombettoni (2001,2003) and corrections to the tunneling and local energy due to the
! presence of the trapping potential

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

!  external k_j
!MODULE K_J_MOD
!  INTERFACE
!     FUNCTION K_J(SITE,BRA,KET)
!       !!       IMPLICIT NONE
!       INTEGER, INTENT(IN)    :: SITE,BRA,KET
!       COMPLEX*16 :: K_J
!     END FUNCTION K_J
!  END INTERFACE
!END MODULE K_J_MOD



PROGRAM TIGH_BINDING
  USE PRINT_MATRIX

  IMPLICIT NONE
  INTEGER INFO, N_BODIES,N_SITES,N_T,D_BARE,I_,J_
  DOUBLE PRECISION TIME,D_H,DT,NORM,MOMENTUM,K_0,TOUT

  EXTERNAL ZVODE_CALL, SCHRODINGER_RHS
  
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: H
  COMPLEX*16, DIMENSION(:),   ALLOCATABLE :: PSI
  DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: E_BARE
  complex*16, dimension(32) :: test
  INTEGER MF

  INFO     = 0

  N_BODIES = 1                          ! number of particles
  N_SITES  = 128                       ! number of lattice sites 
  D_BARE   = D_H(N_SITES,N_BODIES,'B')  ! D_H EVALUATES THE NUMBER OF STATES  

  ALLOCATE(H(D_BARE,D_BARE))
  ALLOCATE(E_BARE(D_BARE))
  ALLOCATE(PSI(D_BARE))

  H = DCMPLX(0.0,0.0)

  N_T  = 256     ! number of time steps
  DT   = 30.0/N_T ! time step
  TIME = 0.0     ! initial time

  
  DO i_=0,N_T
     TIME = i_*DT
     PSI  = 0.0
     CALL TIGHT_BINDING_HAMILTONIAN(SIZE(H,1),TIME,PSI,H,INFO)   
     CALL LAPACK_FULLEIGENVALUES(H,D_BARE,E_BARE,INFO)
     !WRITE(*,*) "#LIST OF INSTANTANEOUS NON-INTERACTING GROUND STATE"
     !WRITE(*,*) "#TIME  POSITION REAL(PSI) IMAG(PSI)"
     DO J_=1,SIZE(E_BARE,1)
        !CALL GET_MOMENTUM(SIZE(E_BARE,1),H(J_,:),MOMENTUM,INFO)
        WRITE(*,*) time,J_,REAL(H(J_,1)),AIMAG(H(J_,1))
     END DO
     WRITE(*,*)
  END DO

  WRITE(*,*)
  WRITE(*,*)

  ! EVALUATE THE SPECTRUM OF THE TIGHT BINDING MODLE
  ! TO SET THE INITIAL CONDITION
  TIME = 0.0
  PSI  = 0.0
  CALL TIGHT_BINDING_HAMILTONIAN(SIZE(H,1),TIME,PSI,H,INFO)  
  CALL LAPACK_FULLEIGENVALUES(H,D_BARE,E_BARE,INFO)
  PSI      = (H(:,1) )!+ H(:,2) + H(:,3))/sqrt(3.0)! INITIAL CONDITION: THE WAVEFUNCTION OF THE GROUND STATE
  NORM     = DOT_PRODUCT(PSI, PSI)
  PSI      = PSI/SQRT(NORM)
  DO J_=1,SIZE(E_BARE,1)
     WRITE(*,*) 0.0,J_,REAL(PSI(J_)),AIMAG(PSI(J_))
  END DO
  WRITE(*,*)
  !WRITE(*,*)
  write(*,*) "# WAVEFUNCTION AS A FUNCTION OF TIME"
  WRITE(*,*) "# TIME  POSITION REAL(PSI) IMAG(PSI)"

  TIME = 0.0
  MF   = 10 ! type of integrator
  INFO = 0
  TOUT = DT

  DO i_=1,N_T
     CALL RKFOURTHORDER(TIGHT_BINDING_HAMILTONIAN,SIZE(E_BARE,1),PSI,TOUT,DT,INFO)
     !CALL RKFOURTHORDER(SIZE(E_BARE,1),PSI,TOUT,DT,INFO)
     !CALL ZVODE_CALL(SCHRODINGER_RHS,SIZE(PSI,1),PSI,TIME,TOUT,MF,INFO)
     
     TOUT = TOUT + DT
     !     IF(MOD(I_,4).EQ.0) THEN
     DO J_=1,SIZE(E_BARE,1)
        WRITE(*,*) TOUT,J_,REAL(PSI(J_)),AIMAG(PSI(J_))
     END DO
     WRITE(*,*)
     !     END IF     
  END DO
   
END PROGRAM TIGH_BINDING


SUBROUTINE SCHRODINGER_RHS(NEQ,TIME,PSI,PSIDOT,RPAR,IPAR)
!      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
  IMPLICIT NONE
  INTEGER,                     INTENT(IN)   :: NEQ
  DOUBLE PRECISION,            INTENT(IN)   :: TIME
  COMPLEX*16, DIMENSION(NEQ),  INTENT(IN)  :: PSI
  COMPLEX*16, DIMENSION(NEQ), INTENT(OUT) :: PSIDOT
  COMPLEX*16,                  INTENT(IN)  :: RPAR
  INTEGER,                     INTENT(IN)  :: IPAR

  EXTERNAL TIGHT_BINDING_HAMILTONIAN
  !YDOT(1) = -RPAR*Y(1)*Y(1)*Y(2)
  !YDOT(2) = RPAR*Y(2)
  COMPLEX*16, DIMENSION(NEQ,NEQ) :: H

  INTEGER INFO
  INFO = 0
  
  CALL TIGHT_BINDING_HAMILTONIAN(NEQ,TIME,PSI,H,INFO)   
  PSIDOT = DCMPLX(0.0,-1.0)*MATMUL(H,PSI)  

END SUBROUTINE SCHRODINGER_RHS
  

