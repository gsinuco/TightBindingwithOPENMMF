

SUBROUTINE Inhomogeneous_Tunneling_B(D_BARE,N_SITES,N_BODIES,STATE,H_J,INFO)
   
  USE CREATIONDESTRUCTION  
  IMPLICIT NONE    
  INTEGER,                               INTENT(IN)    :: D_BARE,N_SITES,N_BODIES
  INTEGER,    DIMENSION(D_BARE,N_SITES), INTENT(IN)    :: STATE
  COMPLEX*16, DIMENSION(D_BARE,D_BARE),  INTENT(OUT)   :: H_J
  INTEGER,                               INTENT(INOUT) :: INFO

  INTEGER, DIMENSION(N_SITES) :: NEW_STATE,STATE_J,STATE_I
    
  DOUBLE PRECISION :: EPSILON_J,DELTA_J,J_REN,K_J

  INTEGER :: N,I_,J_,k_ ,SITE
  N = D_BARE

    
!!$!  write(*,*) D_bare,N_sites
!!$  H_J = 0
!!$  J_=1
!!$  K_=1
!!$
!!$  DO J_=2,D_BARE-1
!!$     SITE = k_ - (N_SITES-1)/2
!!$     CALL TRAP_EFFECTS(SITE,EPSILON_J,DELTA_J,J_REN,K_J,INFO)
!!$     H_J(J_,J_+1) = DELTA_J
!!$     H_J(J_,J_-1) = DELTA_J
!!$  END DO
!!$  CALL TRAP_EFFECTS(1,EPSILON_J,DELTA_J,J_REN,K_J,INFO)
!!$  H_J(1,2)      = DELTA_J
!!$  H_J(1,D_BARE) = DELTA_J
!!$  CALL TRAP_EFFECTS(D_BARE,EPSILON_J,DELTA_J,J_REN,K_J,INFO)
!!$  H_J(D_BARE,D_BARE-1) = DELTA_J
!!$  H_J(D_BARE,1)        = DELTA_J
  

  DO k_=1,N_SITES-1 ! loop though all sites
     SITE = k_ - (N_SITES-1)/2
     CALL TRAP_EFFECTS(SITE,EPSILON_J,DELTA_J,J_REN,K_J,INFO)
     DO J_=1,N      ! loop through all states
        DO I_=J_+1,N ! and evaluate the forward tunneling matrix
           !WRITE(*,*) K_,J_,I_,STATE(J_,:) 
            !WRITE(*,*) STATE(I_,:) 
           STATE_J = STATE(J_,:) 
            STATE_I = STATE(I_,:) 
            NEW_STATE = TUNNELING_(k_,STATE_J)
            !write(*,*) NEW_STATE,dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I))
            IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
               !WRITE(*,*) K_J(k_,J_,I_)
                H_J(I_,J_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))*DELTA_J
                H_J(J_,I_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))*DELTA_J
            END IF           
        END DO
     END DO

     DO J_=1,N
         STATE_J = STATE(J_,:) 
         STATE_I = STATE(J_,:) 
         NEW_STATE = TUNNELING_(k_,STATE_J)
         IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
            H_J(J_,J_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))*DELTA_J
         END IF           
     END DO
  END DO

  !PERIODIC BOUNDARY CONDITIONS
  H_J(N_SITES,1) = H_J(1,2)
  H_J(1,N_SITES) = H_J(1,2)

END SUBROUTINE Inhomogeneous_Tunneling_B

SUBROUTINE Inhomogeneous_EnergyShift_B(D_BARE,N_SITES,N_BODIES,STATE,H_J,INFO)
   
  USE CREATIONDESTRUCTION  
  IMPLICIT NONE    
  INTEGER,                               INTENT(IN)    :: D_BARE,N_SITES,N_BODIES
  INTEGER,    DIMENSION(D_BARE,N_SITES), INTENT(IN)    :: STATE
  COMPLEX*16, DIMENSION(D_BARE,D_BARE),  INTENT(OUT)   :: H_J
  INTEGER,                               INTENT(INOUT) :: INFO

  INTEGER, DIMENSION(N_SITES) :: NEW_STATE,STATE_J,STATE_I
    
  DOUBLE PRECISION :: EPSILON_J,DELTA_J,J_REN,K_J

  INTEGER :: N,I_,J_,k_,site 
  N = D_BARE

    
!  write(*,*) D_bare,N_sites
  H_J = 0
  J_=1
  K_=1
  
  DO k_=1,N_SITES-1 ! loop though all sites     
     SITE = k_ - (N_SITES-1)/2
     CALL TRAP_EFFECTS(SITE,EPSILON_J,DELTA_J,J_REN,K_J,INFO)
     H_J(k_,k_) = EPSILON_J
  END DO

  !PERIODIC BOUNDARY CONDITIONS
  !H_J(N_SITES,1) = H_J(N_SITES,N_SITES)
  !H_J(1,N_SITES) = H_J(1,1)

END SUBROUTINE Inhomogeneous_EnergyShift_B



SUBROUTINE Tunneling_B(D_BARE,N_SITES,N_BODIES,STATE,H_J,INFO)
   
  USE CREATIONDESTRUCTION
  IMPLICIT NONE    
  INTEGER,                               INTENT(IN)    :: D_BARE,N_SITES,N_BODIES
  INTEGER,    DIMENSION(D_BARE,N_SITES), INTENT(IN)    :: STATE
  COMPLEX*16, DIMENSION(D_BARE,D_BARE),  INTENT(OUT)   :: H_J
  INTEGER,                               INTENT(INOUT) :: INFO

  INTEGER, DIMENSION(N_SITES) :: NEW_STATE,STATE_J,STATE_I
    
  INTEGER :: N,I_,J_,k_ 
  N = D_BARE

!  write(*,*) D_bare,N_sites
  H_J = 0
  J_=1
  K_=1
  DO k_=1,N_SITES-1 ! loop though all sites
     DO J_=1,N      ! loop through all states
        DO I_=J_+1,N ! and evaluate the forward tunneling matrix
           !WRITE(*,*) K_,J_,I_,STATE(J_,:) 
            !WRITE(*,*) STATE(I_,:) 
            STATE_J = STATE(J_,:) 
            STATE_I = STATE(I_,:) 
            NEW_STATE = TUNNELING_(k_,STATE_J)
!            write(*,*) NEW_STATE,dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I))
            IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
                H_J(I_,J_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))
                H_J(J_,I_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))
            END IF           
        END DO
     END DO

     DO J_=1,N
         STATE_J = STATE(J_,:) 
         STATE_I = STATE(J_,:) 
         NEW_STATE = TUNNELING_(k_,STATE_J)
         IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
            H_J(J_,J_) = SQRT(1.0*STATE_J(k_)*(STATE_J(k_+1)+1))
         END IF           
     END DO
  END DO

  !PERIODIC BOUNDARY CONDITIONS
  H_J(N_SITES,1) = H_J(N_SITES,N_SITES-1)
  H_J(1,N_SITES) = H_J(1,2)

END SUBROUTINE Tunneling_B

SUBROUTINE Onsite_meanfield_B(D_BARE,PSI,H_U,INFO)

    USE CREATIONDESTRUCTION
  
    IMPLICIT NONE    
    INTEGER,                               INTENT(IN)    :: D_BARE
    COMPLEX*16, DIMENSION(D_BARE),         INTENT(IN)    :: PSI
    COMPLEX*16, DIMENSION(D_BARE,D_BARE),  INTENT(OUT)   :: H_U
    INTEGER,                               INTENT(INOUT) :: INFO

    INTEGER :: N,I_,J_,k_ 
    N = D_BARE

    H_U = 0
    DO J_=1,N
       H_U(J_,J_) = ABS(PSI(J_))**2
    END DO

END SUBROUTINE Onsite_meanfield_B

SUBROUTINE Onsite_twobody_B(D_BARE,N_SITES,N_BODIES,STATE,H_U,INFO)

    USE CREATIONDESTRUCTION
  
    IMPLICIT NONE    
    INTEGER,                               INTENT(IN)    :: D_BARE,N_SITES,N_BODIES
    INTEGER,    DIMENSION(D_BARE,N_SITES), INTENT(IN)    :: STATE
    COMPLEX*16, DIMENSION(D_BARE,D_BARE),  INTENT(OUT)   :: H_U
    INTEGER,                               INTENT(INOUT) :: INFO

    DOUBLE PRECISION, DIMENSION(N_SITES) :: NEW_STATE
    
    
    INTEGER :: N,I_,J_,k_ 
    N = D_BARE

    H_U = 0
    DO k_=1,N_SITES ! loop though all sites
        DO J_=1,N
            H_U(J_,J_) = STATE(J_,k_)
            H_U(J_,J_) = H_U(J_,J_)*(H_U(J_,J_)-1.0)
        END DO
    END DO

END SUBROUTINE Onsite_twobody_B

SUBROUTINE Tunneling_F(D_BARE,N_SITES,N_BODIES,STATE,H_J,INFO)
   
  USE CREATIONDESTRUCTION
  IMPLICIT NONE    
  INTEGER,                               INTENT(IN)    :: N_SITES
  INTEGER,    DIMENSION(2),              INTENT(IN)    :: D_BARE,N_BODIES
  INTEGER,    DIMENSION(D_BARE(1)+D_BARE(2),N_SITES),      INTENT(IN)    :: STATE
  COMPLEX*16, DIMENSION(D_BARE(1)+D_BARE(2),D_BARE(1)+D_BARE(2)), INTENT(OUT)   :: H_J
  INTEGER,                               INTENT(INOUT) :: INFO

  INTEGER, DIMENSION(N_SITES) :: NEW_STATE,STATE_J,STATE_I
    
  COMPLEX*16, DIMENSION(D_BARE(1),D_BARE(1)) :: T_UP
  COMPLEX*16, DIMENSION(D_BARE(2),D_BARE(2)) :: T_DOWN
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: T
  
  
  INTEGER :: N,I_,J_,k_,l_
  INTEGER, DIMENSION(3) :: N_
  
  N_(1) = D_BARE(1)
  N_(2) = D_BARE(2)
  N_(3) = D_BARE(1)*D_BARE(2)
  ALLOCATE(T(N_(3),N_(3)))
  
  T_UP   = 0
  T_DOWN = 0
  H_J    = 0
  T      = 0
  DO k_=1,N_SITES ! loop though all sites
    DO l_=1,2 ! loop through spin up and spin down
        N = D_BARE(l_)
        DO J_=1,N ! Nested loop through all spin up/down states
            DO I_=J_,N
                STATE_J = STATE(J_ + (l_-1)*D_BARE(1),:) 
                STATE_I = STATE(I_ + (l_-1)*D_BARE(1),:) 
                NEW_STATE = TUNNELING_F_(k_,STATE_J)
                
                
                
                
                !write(*,*) NEW_STATE,dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I))
                IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
                    IF(l_.EQ.1) THEN
                        T_UP(I_,J_) = 1.0
                        T_UP(J_,I_) = 1.0
                    ELSE
                        T_DOWN(I_,J_) = 1.0
                        T_DOWN(J_,I_) = 1.0
                    END IF  

                END IF           
            END DO
        END DO

        DO J_=1,N
            STATE_J = STATE(J_+(l_-1)*D_BARE(1),:) 
            STATE_I = STATE(J_+(l_-1)*D_BARE(1),:) 
            NEW_STATE = TUNNELING_F_(k_,STATE_J)
            IF(dot_product((NEW_STATE-STATE_I),(NEW_STATE-STATE_I)).EQ.0) THEN
                IF(l_.EQ.1) THEN
                    T_UP(J_,J_) = 1.0
                ELSE
                    T_DOWN(J_,J_) = 1.0
                END IF
            END IF           
        END DO
    END DO  
    CALL TENSORMULT(N_,T_UP,T_DOWN,T,INFO)
    H_J = H_J + T
  END DO

END SUBROUTINE Tunneling_F


SUBROUTINE Onsite_twobody_F(D_BARE,N_SITES,N_BODIES,STATE,H_U,INFO)
   
  USE CREATIONDESTRUCTION
  IMPLICIT NONE    
  INTEGER,                               INTENT(IN)    :: N_SITES
  INTEGER,    DIMENSION(2),              INTENT(IN)    :: D_BARE,N_BODIES
  INTEGER,    DIMENSION(D_BARE(1)+D_BARE(2),N_SITES),      INTENT(IN)    :: STATE
  COMPLEX*16, DIMENSION(D_BARE(1)+D_BARE(2),D_BARE(1)+D_BARE(2)), INTENT(OUT)   :: H_U
  INTEGER,                               INTENT(INOUT) :: INFO

  INTEGER, DIMENSION(N_SITES) :: NEW_STATE,STATE_J,STATE_I    

  
  INTEGER :: N,I_,J_,k_,l_
  
  COMPLEX*16, DIMENSION(D_BARE(1),D_BARE(1)) :: T_UP
  COMPLEX*16, DIMENSION(D_BARE(2),D_BARE(2)) :: T_DOWN
  COMPLEX*16, DIMENSION(:,:), ALLOCATABLE :: T

  INTEGER, DIMENSION(3) :: N_
  
  H_U = 0
  T_UP   = 0
  T_DOWN = 0
  N_(1) = D_BARE(1)
  N_(2) = D_BARE(2)
  N_(3) = D_BARE(1)*D_BARE(2)
  ALLOCATE(T(N_(3),N_(3)))
    
  !write(*,*) N_ 
  T = 0.0
  DO l_=1,N_(2)
      DO k_=1,N_(1)
          !write(*,*) N_,(l_-1)*N_(1)+k_,(l_-1)*N_(1)+k_,DOT_PRODUCT(STATE(k_,:),STATE(N_(1)+l_,:))
        T((l_-1)*N_(1)+k_,(l_-1)*N_(1)+k_) = DOT_PRODUCT(STATE(k_,:),STATE(N_(1)+l_,:))
    END DO
  END DO
  H_U = T
END SUBROUTINE Onsite_twobody_F


SUBROUTINE GET_MOMENTUM(N,PSI,MOMENTUM,INFO)
  !# The energy eigenstates are also eigenstates of the translation vector T_n, 
  !# where T_n psi_k(x) = psi_k(a*n+x). With T_n psi_k(x) = exp(ikna)psi_k(x) 
  !# => k = arg(pis_k(x)/psi(k+a))/a. Here, a = 1, f(a) = f[m]
  !def get_k(v_i, m):
  !    return np.angle(v_i[0]/ v_i[m])
  
  IMPLICIT NONE
  INTEGER,                         INTENT(IN)    :: N
  INTEGER,                         INTENT(INOUT) :: INFO
  COMPLEX*16,       DIMENSION(N,N),INTENT(IN)    :: PSI
  DOUBLE PRECISION, DIMENSION(N),  INTENT(OUT)   :: MOMENTUM

  INTEGER i
  COMPLEX*16 PSI_RATIO
  
  IF(INFO.EQ.0) THEN
     DO i=1,SIZE(MOMENTUM,1)
        PSI_RATIO   = PSI(1,i)/psi(2,i)
        MOMENTUM(i) = ATAN(AIMAG(PSI_RATIO)/REAL(PSI_RATIO))
     END DO
  ELSE
     WRITE(*,*) "INFO != 0"   
     WRITE(*,*) "WE COULD NOT EVALUATE THE MOMENTUM"
  END IF


END SUBROUTINE GET_MOMENTUM

