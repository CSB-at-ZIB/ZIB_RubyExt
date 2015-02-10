      PROGRAM TNLSCO
      IMPLICIT DOUBLEPRECISION(S)
C
C     ____________________________________________________________
C
C     Testexample for Nlscon: Least Squares Approximation of
C     measurement points with an Exponential of Negative Square
C     curve.
C
C*  Written by        L. Weimann 
C*  Purpose           Testexample for code NLSCON
C*  Version           2.3
C*  Revision          January 1992
C*  Latest Change     June 1992
C*  Library           CodeLib
C*  Code              Fortran 77, Double Precision
C*  Environment       Standard Fortran 77 environment on PC's,
C                     workstations and hosts.
C
C     ____________________________________________________________
C
      INTEGER IRW
      PARAMETER (IRW=400)
      INTEGER IIW
      PARAMETER (IIW=60)
      INTEGER N,M,MFIT,I
      DOUBLE PRECISION EPS,H
      INTEGER IOPT(50)
      INTEGER IERR,IFAIL
      DOUBLE PRECISION X(3),XUSCAL(3),FOBS(11),FWEIGH(11),RW(IRW)
      INTEGER IW(IIW)
      REAL STIME, ETIME, CPTIME
      EXTERNAL F
      EXTERNAL DF
C:    Begin
      OPEN(2,FILE='nlscon.dat')
      OPEN(9,FILE='nlscon.out')
      WRITE(6,*) ' monitor: nlscon.out , data: nlscon.dat '
C     Number of parameters to be estimated
      N = 3
C     Number of observations ( data points fobs ) + number of
C       equality constraints
      M = 12
C     Number of observations
      MFIT = 11
C     Exact parameters of the model function for the computation
C       of the data points FOBS(I)( data to be estimated in real
C       life problems )
      X(1)=100.0D0
      X(2)=10.0D0
      X(3)=10.0D0
C     Evaluation of the model function for these parameters ( data
C       gained by measurements in real life problems - measurement
C       errors ( less or equal 0.01 ) are simulated by sin function
C       term )
80    FORMAT(/,' Simulated  experimential  data  :')
      WRITE(9,80)
      DO 81 I=1,M-1
        H =(2.0D0*(DBLE(I)-1.0D0)-X(3))/X(2)
        FOBS(I)=X(1)*DEXP(-H*H/2.0D0)+DSIN(10.0D0*I)*1.0D-2
        WRITE(9,*)FOBS(I)
81    CONTINUE
82    FORMAT(//)
      WRITE(9,82)
C
      EPS = 1.0D-4
      DO 710 I=1,50
        IOPT(I)=0
710    CONTINUE
      DO 711 I=1,IIW
        IW(I)=0
711   CONTINUE
      DO 712 I=1,IRW
        RW(I)=0.0D0
712   CONTINUE
C     Execution mode: 0=Standard Mode, 1=Stepwise mode
      IOPT(2)=1
C       Jacobian: 0=(same as value 3)
C                 1=supplied by user routine JAC
C                 2=computed by numerical differentation (no feedback) 
C                 3=computed by numerical differentation (with feedback)
      IOPT(3)=1
C     A posteriori statistical analysis: 0 = no, 1 = yes
      IOPT(21)=1
C     Broyden updates: 0 = inhibit, 1=allow
      IOPT(32)=0
C     Problem classification:
C     1 = linear , 2 = mildly nonlinear  3 = highly nonlinear
      IOPT(31)=3
C     Automatic row scaling of linear system (constraints)
C     and user scaling (measurements):
C     0 = allowed , 1 = inhibited
      IOPT(35)=0
C       Set MPRERR, MPRMON, MPRSOL, MPRTIM, MPRSTA
        IOPT(11)=3
        IOPT(13)=3
        IOPT(15)=2
        IOPT(19)=1
        IOPT(22)=1
C       Set print units LUERR, LUMON, LUSOL, LUTIM
        IOPT(12)=9
        IOPT(14)=9
        IOPT(16)=2
        IOPT(20)=9
C     Solution output format:
C     0=standard format, 1= GRAZIL readable output
      IOPT(46)=0
C     Override maximum allowed number of iterations:
      IW(31)=200
C     Override initial pseudo-rank:
C     IW(32)=N
C     Override starting damping factor:
C     RW(21)=1.0D0
C     Override minimal allowed damping factor:
C     RW(22)=1.0D-3
C     Override rank1-decision parameter SIGMA:
C     RW(23)=2.0D0
C     Override maximum permitted subcondition for DECCON:
      RW(25)= 1.0D+16
C     Initial guess of parameters to be estimated
      X(1)=1.0D0
      X(2)=2.0D0
      X(3)=5.0D0
      DO 15 I=1,N
        XUSCAL(I) = 0.0
15    CONTINUE
      DO 20 I=1,MFIT
        FWEIGH(I)=0.0
20    CONTINUE
      IERR=-1
      I=0
      CALL ZIBSEC(STIME,IFAIL)
31    IF (IERR.EQ.-1) THEN
        CALL NLSCON(N,M,MFIT,F,DF,X,XUSCAL,FOBS,FWEIGH,EPS,IOPT,IERR,
     $  IIW,IW,IRW,RW)
C       Clear workspace declared not to be used
        NIFREE=IW(16)
        DO 311 K=NIFREE,IIW
          IW(K)=0
311     CONTINUE
        NRFREE=IW(17)
        DO 312 K=NRFREE,IRW
          RW(K)=0.0D0
312     CONTINUE
        I=I+1
32       FORMAT(' Returned from call ',I4,' of NLSCON')
        WRITE(9,32)I
C       IOPT(2)=0
        GOTO 31
      ENDIF
      CALL ZIBSEC(ETIME,IFAIL)
      CPTIME = ETIME-STIME
73    FORMAT(//,1X,'Time ','used ','=',F9.3,1X,'Sec',//,66('*'),
     $  /)
      WRITE(9,73)CPTIME
      END
      SUBROUTINE F(N,M,MCON,X,FX,IFLAG)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N,M,MCON
      DOUBLE PRECISION X(N),FX(M)
C:    End Parameter
      INTEGER I
      DOUBLE PRECISION H
C:    Begin
      FX(1)=X(1)-X(2)**2-X(3)**2+100.0D0
      DO 84 I=2,M
        H =(2.0D0*(DBLE(I-1)-1.0D0)-X(3))/X(2)
        FX(I)=X(1)*DEXP(-H*H/2.0D0)
84    CONTINUE
      RETURN
      END
      SUBROUTINE DF(N,M,MCON,X,DFX,IFLAG)
      IMPLICIT DOUBLEPRECISION(S)
      INTEGER N,M,MCON
      DOUBLE PRECISION X(N),DFX(M,N)
C:    End Parameter
      INTEGER I
      DOUBLE PRECISION H,EXPH,HEXPH
C:    Begin
      DFX(1,1)=1.0D0
      DFX(1,2)=-2.0D0*X(2)
      DFX(1,3)=-2.0D0*X(3)
      DO 85 I=2,M
        H =(2.0D0*(DBLE(I-1)-1.0D0)-X(3))/X(2)
        EXPH = DEXP(-H*H/2)
        HEXPH = X(1)*EXPH*H/X(2)
        DFX(I,1)=EXPH
        DFX(I,2)=HEXPH*H
        DFX(I,3)=HEXPH
85    CONTINUE
      RETURN
      END
