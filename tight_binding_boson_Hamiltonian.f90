
SUBROUTINE TIGHT_BINDING_HAMILTONIAN(N,time,H,INFO)
  !USE ATOMIC_PROPERTIES
  !USE TYPES
  USE SUBINTERFACE
  USE SUBINTERFACE_LAPACK
  !USE FLOQUETINITINTERFACE
  !USE ARRAYS 


  IMPLICIT NONE
  INTEGER,                    INTENT(IN)    :: N
  DOUBLE PRECISION,           INTENT(IN)    :: time
  COMPLEX*16, DIMENSION(N,N), INTENT(OUT)   :: H
  INTEGER,                    INTENT(INOUT) :: INFO
  

  !TYPE(MODE),       DIMENSION(:),   ALLOCATABLE :: FIELDS
  !TYPE(ATOM)                                       ID
  INTEGER,          DIMENSION(:),   ALLOCATABLE :: MODES_NUM
  INTEGER                                          TOTAL_FREQUENCIES,D_BARE
  INTEGER                                          m,INDEX0,r
  DOUBLE PRECISION, DIMENSION(:),   ALLOCATABLE :: ENERGY,E_FLOQUET,E_BARE,MOMENTUM
  COMPLEX*16,       DIMENSION(:),   ALLOCATABLE :: PSI_K
  COMPLEX*16,       DIMENSION(:,:), ALLOCATABLE :: H_J,H_U,U_F,H_BARE
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: P_AVG
  INTEGER         , DIMENSION(:,:), ALLOCATABLE :: states_occ
  DOUBLE PRECISION                              :: T1,T2

  INTEGER index,N_BODIES,N_SITES,i_,M_
  DOUBLE PRECISION t,u,t_driv,omega,D_H

  !OPEN(UNIT=3,FILE="ManybodyHubbard_Bosons.dat",ACTION="WRITE")
  
  IF (INFO .EQ. 0) THEN
     N_BODIES = 1
     
     N_SITES  = N
     
     
     
     D_BARE = D_H(N_SITES,N_BODIES,'B')  ! D_H EVALUATES THE NUMBER OF STATES  
     write(*,*) '# Number of lattice sites:         ', N_SITES
     write(*,*) '# Number of particles:             ', N_BODIES
     write(*,*) '# Number of states Bosonic states  ', D_BARE
     write(*,*)
     write(*,*)
     
     
         
     ALLOCATE(E_BARE(D_BARE))             ! STORE THE ENERGY SPECTRUM
     ALLOCATE(H_BARE(D_BARE,D_BARE))      ! STORE THE BARE HAMILTONIAN
     ALLOCATE(H_J(D_BARE,D_BARE))         ! STORE THE TUNNENING MATRIX
     ALLOCATE(H_U(D_BARE,D_BARE))         ! STORES THE ONSITE INTERACTION
     ALLOCATE(MOMENTUM(D_BARE))           ! STORES THE MOMENTUM ASSOCIATED WITH AN EIGENENERGY
     ALLOCATE(PSI_K(D_BARE))              ! STORES THE MOMENTUM ASSOCIATED WITH AN EIGENENERGY
     ALLOCATE(MODES_NUM(2))               ! NUMBER OF DRIVING MODES
     
     
     ALLOCATE(states_occ(D_BARE,N_SITES)) ! STORES THE STATES USING OCCUPATION NUMBER 
     
     ! CREATE THE BASIS OF STATES
     CALL Manybody_basis(D_BARE,N_SITES,N_BODIES,'B',states_occ,INFO)
     !call write_matrix_int(states_occ)
     
     
     
     !  ! EVALUATE THE TUNNELING TERM OF THE HAMILTONIAN
     
     
     
     
     
     
     
     
     
     CALL Tunneling_B(D_BARE,N_SITES,N_BODIES,STATES_OCC,H_J,INFO)
     !  CALL WRITE_MATRIX(abs(h_j))
     !  !! EVALUATE THE ON-SITE INTERACTION
     
     CALL Onsite_twobody_B(D_BARE,N_SITES,N_BODIES,STATES_OCC,H_U,INFO)
     !  !CALL WRITE_MATRIX(abs(h_u))
     !
     ! HUBBARD MODEL PARAMETERS
     t      = 0.1 - 0.2*TIME**2
     u      = 0.0
     
     write(*,*) '# Hubbard parameters (t,u):',t,u
     !IF(TIME.EQ.0) H =  -t*H_J + u*H_U     
     !IF(TIME.GT.0) H =   t*H_J + u*H_U     
     H =  -t*H_J + u*H_U     
     
     ! EVALUATE THE SPECTRUM OF THE STATIC HAMILTONIAN
     !H_BARE = FIELDS(1)%V
     !CALL LAPACK_FULLEIGENVALUES(H_BARE,D_BARE,E_BARE,INFO)
     !CALL WRITE_MATRIX(real(H_BARE))
     !CALL WRITE_MATRIX(aimag(H_BARE))
     
     !The eigenvectors are the columns of H_BARE, i.e.
     !<n|e_i> = H_BARE(n,i)
     !where |n> is the n-th lattice site. This is because:
     ! |e_i> = SUM_n H_BARE(n,i) |n>
     !CALL GET_MOMENTUM(SIZE(H_BARE,1),H_BARE,MOMENTUM,INFO)
     !CALL FFTW_1D(SIZE(H_BARE,1),H_BARE(:,64),PSI_K,"FORWARD",INFO)
     !DO I_=1,SIZE(E_BARE)
     !MOMENTUM(I_) = (-1)**(I_-1)
     !write(*,*) I_,E_BARE(I_),abs(H_BARE(I_,2)),AIMAG(H_BARE(I_,2)),REAL(H_BARE(I_,3)),AIMAG(H_BARE(I_,3)),&
     !     &     real(PSI_K(I_))
     
     !END DO
     !WRITE(*,*)

  ELSE
     WRITE(*,*) "INFO != 0"
     write(*,*) "nothing is done!"
  END IF

END SUBROUTINE TIGHT_BINDING_HAMILTONIAN
  

