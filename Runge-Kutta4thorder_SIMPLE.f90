SUBROUTINE RKFOURTHORDER(N,y,t,dt,INFO)
  
  IMPLICIT NONE
  INTEGER,                          INTENT(IN)    :: N
  COMPLEX*16,       DIMENSION(N),   INTENT(INOUT) :: y
  DOUBLE PRECISION,                 INTENT(INOUT) :: t
  DOUBLE PRECISION,                 INTENT(IN)    :: dt
  INTEGER,                          INTENT(INOUT) :: INFO
  !DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: beta_ijkl
  !DOUBLE PRECISION, DIMENSION(:), INTENT(IN)    :: beta_ijkl_map
  
  COMPLEX*16, DIMENSION(:), ALLOCATABLE :: K1,K2,K3,K4,y_aux
  COMPLEX*16, DIMENSION(N,N) :: H
  
  
  INTEGER i,dim
  
  IF (INFO.EQ.0) THEN
     ALLOCATE(K1(SIZE(Y)))
     ALLOCATE(K2(SIZE(Y)))
     ALLOCATE(K3(SIZE(Y)))
     ALLOCATE(K4(SIZE(Y)))
     ALLOCATE(y_aux(SIZE(Y)))
     
     !write(*,*) "k1"
     !K1 = f(t,y,beta_ijkl,beta_ijkl_map) 
     CALL TIGHT_BINDING_HAMILTONIAN(N,t,Y,H,INFO)
     y_aux = -DCMPLX(0.0,1.0)*MATMUL(H,Y)
     K1 = Y_AUX !f(t,y)
     y_aux = y +K1*dt*0.5
     
     !write(*,*) "k2"
     !K2 = f(t+0.5*dt,y_aux,beta_ijkl,beta_ijkl_map)
     CALL TIGHT_BINDING_HAMILTONIAN(N,t+0.5*dt,y_aux,H,INFO)
     y_aux = -DCMPLX(0.0,1.0)*MATMUL(H,y_aux)
     K2 = y_aux!f(t+0.5*dt,y_aux)
     y_aux = y +K2*dt*0.5
     
     !write(*,*) "k3"
     !K3 = f(t+0.5*dt,y_aux,beta_ijkl,beta_ijkl_map)
     CALL TIGHT_BINDING_HAMILTONIAN(N,t+0.5*dt,Y_aux,H,INFO)
     y_aux = -DCMPLX(0.0,1.0)*MATMUL(H,y_aux)
     K3 = y_aux!f(t+0.5*dt,y_aux)
     y_aux = y +K3*dt
     
     !write(*,*) "k4"
     !K4 = f(t+dt,y_aux,beta_ijkl,beta_ijkl_map)
     CALL TIGHT_BINDING_HAMILTONIAN(N,t+dt,Y_aux,H,INFO)
     y_aux = -DCMPLX(0.0,1.0)*MATMUL(H,y_aux)
     K4 = y_aux!f(t+dt,y_aux)
     y = y + dt*(K1 + 2.0*K2 + 2.0*K3 + K4)/6.0
     
     t = t + dt;
     !write(*,*) t,dt  
  ELSE
     WRITE(*,*) "INFO != 0"
     WRITE(*,*) "WARNING: NOTHING IS DONE"
  END IF
END SUBROUTINE RKFOURTHORDER
  
