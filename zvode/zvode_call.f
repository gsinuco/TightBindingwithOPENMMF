C     EXAMPLE PROBLEM
C     
C     The program below uses ZVODE to solve the following system of 2 ODEs:
C     dw/dt = -i*w*w*z, dz/dt = i*z; w(0) = 1/2.1, z(0) = 1; t = 0 to 2*pi.
C     Solution: w = 1/(z + 1.1), z = exp(it).  As z traces the unit circle,
C     w traces a circle of radius 10/2.1 with center at 11/2.1.
C     For convenience, Main passes RPAR = (imaginary unit i) to FEX and JEX.
C     

c$$$       PROGRAM zvode_test
c$$$       IMPLICIT NONE
c$$$       COMPLEX*16, DIMENSION(:), ALLOCATABLE :: Y
c$$$       INTEGER NEQ
c$$$       DOUBLE PRECISION T,DT,TOUT
c$$$       INTEGER MF,INFO,IOUT
c$$$
c$$$       EXTERNAL FEX
c$$$
c$$$       NEQ = 2
c$$$       ALLOCATE(Y(NEQ))
c$$$
c$$$       T  = 0
c$$$       DT = 0.1570796326794896D0
c$$$       MF = 10
c$$$       INFO =0
c$$$
c$$$       Y(1)  = 1.0D0/2.1D0
c$$$       Y(2)  = 1.0D0
c$$$  
c$$$       TOUT = DT
c$$$       WRITE(6,10)
c$$$ 10    FORMAT('   t',11X,'w',26X,'z')
c$$$       DO 41 IOUT = 1,40
c$$$          CALL ZVODE_CALL(FEX,SIZE(Y),Y,T,TOUT,MF,INFO)
c$$$ 41       TOUT = TOUT + DT     
c$$$
c$$$       STOP
c$$$       END 

       SUBROUTINE ZVODE_CALL(FEX,NEQ,Y,T,DT,MF,INFO)
       IMPLICIT NONE
       INTEGER,                    INTENT(IN)    :: NEQ
       COMPLEX*16, DIMENSION(NEQ), INTENT(INOUT) :: Y
       DOUBLE PRECISION,           INTENT(IN)    :: T,DT
       INTEGER,                    INTENT(IN)    :: MF
       INTEGER,                 INTENT(INOUT)    :: INFO

       EXTERNAL FEX, JEX
       DOUBLE COMPLEX,   DIMENSION(:), ALLOCATABLE :: ZWORK
       DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: RWORK
       INTEGER,          DIMENSION(:), ALLOCATABLE :: IWORK

       DOUBLE COMPLEX   :: RPAR, WTRU, ERR, TOUT
       DOUBLE PRECISION :: ABERR, AEMAX, ATOL, RTOL
       INTEGER          :: IOPT,IPAR,ISTATE,ITASK,ITOL,IOUT
       INTEGER          :: LIW,LRW,LZW,ML,MU
       
       TOUT  = DT
       ITOL  = 1
       RTOL  = 1.D-9
       ATOL  = 1.D-8
       ITASK = 1
       ISTATE = 1
       IOPT  = 0
       RPAR  = DCMPLX(0.0D0,1.0D0)
       AEMAX = 0.0D0

       IF(MF.EQ.10) LZW = 15*NEQ                      !              for MF = 10,
       IF(MF.EQ.21 .OR. MF.EQ.22) LZW = 8*NEQ + 2*NEQ**2            !for MF = 21 or 22,
       IF(MF.EQ.24 .OR. MF.EQ.25) LZW = 10*NEQ + (3*ML + 2*MU)*NEQ  !for MF = 24 or 25.
       ALLOCATE(ZWORK(LZW))

       LRW   = 20 + NEQ
       ALLOCATE(RWORK(LRW))

       IF(MF.EQ.10)                LIW = 30        !for MF = 10,
       IF(MF.GE.21 .AND. MF.LE.25) LIW = 30 + NEQ  !for MF = 21, 22, 24, or 25.
       ALLOCATE(IWORK(LIW))

!         CALL ZVODE(TIGHT_BINDING_HAMILTONIAN,D_BARE,PSI,T,DT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,
!     c        ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
         CALL ZVODE(FEX,NEQ,Y,T,TOUT,ITOL,RTOL,ATOL,ITASK,ISTATE,IOPT,
     c        ZWORK,LZW,RWORK,LRW,IWORK,LIW,JEX,MF,RPAR,IPAR)
         WTRU = 1.0D0/DCMPLX(COS(T) + 1.1D0, SIN(T))
         ERR = Y(1) - WTRU
         ABERR = ABS(DREAL(ERR)) + ABS(DIMAG(ERR))
         AEMAX = MAX(AEMAX,ABERR)
!         WRITE(6,20) T, DREAL(Y(1)),DIMAG(Y(1)), DREAL(Y(2)),DIMAG(Y(2))
! 20      FORMAT(F9.5,2X,2F12.7,3X,2F12.7)
         IF (ISTATE .LT. 0) THEN
            WRITE(6,30) ISTATE
 30         FORMAT(//'***** Error halt.  ISTATE =',I3)
            STOP
         ENDIF
 !        WRITE(6,50) IWORK(11), IWORK(12), IWORK(13), IWORK(20),
 !    c       IWORK(21), IWORK(22), IWORK(23), AEMAX
 !50      FORMAT(/' No. steps =',I4,'   No. f-s =',I5,
 !    c        '   No. J-s =',I4,'   No. LU-s =',I4/
 !    c        ' No. nonlinear iterations =',I4/
 !    c        ' No. nonlinear convergence failures =',I4/
 !    c        ' No. error test failures =',I4/
 !    c        ' Max. abs. error in w =',D10.2)
!         STOP
       END
C     
!      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR)
!      DOUBLE COMPLEX Y(NEQ), YDOT(NEQ), RPAR
!      DOUBLE PRECISION T
!      YDOT(1) = -RPAR*Y(1)*Y(1)*Y(2)
!      YDOT(2) = RPAR*Y(2)
!      RETURN
!      END
C     
      SUBROUTINE JEX (NEQ, T, Y, ML, MU, PD, NRPD, RPAR, IPAR)
      DOUBLE COMPLEX Y(NEQ), PD(NRPD,NEQ), RPAR
      DOUBLE PRECISION T
      PD(1,1) = -2.0D0*RPAR*Y(1)*Y(2)
      PD(1,2) = -RPAR*Y(1)*Y(1)
      PD(2,2) = RPAR
      RETURN
      END
