/*
C
C*    External subroutines (to be supplied by the user)
C     =================================================
C 
C     (Caution: Arguments declared as (input) must not
C               be altered by the user subroutines ! )
C
C     FCN(N,M,MCON,X,F,IFAIL) 
C                      Ext    Function subroutine
C       N              Int    Number of vector components (input)
C       M              Int    Number of measurement-vector
C                             components plus number of equality
C                             constraints (input)
C       MCON           Int    Number of equality constraints (input)
C       X(N)           Dble   Vector of parameters (input)
C       F(M)           Dble   Vector of equality constraints and
C                             measurement fitting values -
C                             the first MCON-components belonging
C                             to the constraints (output).
C       IFAIL          Int    FCN evaluation-failure indicator. (output)
C                             On input:  Has always value 0 (zero).
C                             On output: Indicates failure of FCN eval-
C                                uation, if having a value <= 2.
C                             If <0: NLSCON will be terminated with 
C                                    error code = 82, and IFAIL stored
C                                    to IWK(23).
C                             If =1: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced to it's half.
C                             If =2: A new trial Newton iterate will
C                                    computed, with the damping factor
C                                    reduced by a reduct. factor, which 
C                                    must be output through F(1) by FCN,
C                                    and it's value must be >0 and < 1.
C                             Note, that if IFAIL = 1 or 2, additional
C                             conditions concerning the damping factor,
C                             e.g. the minimum damping factor or the
C                             bounded damping strategy may also influ-
C                             ence the value of the reduced damping 
C                             factor.
C
C
C     JAC(N,M,MCON,X,DFDX,IFAIL)
C                        Ext    Jacobian matrix subroutine
C       N                  Int    See parameter N of FCN above (input)
C       M                  Int    See parameter M of FCN above (input)
C       MCON               Int    See parameter MCON of FCN above 
C                                 (input)
C       X(N)               Dble   See parameter X of FCN above (input)
C       DFDX(M,N)          Dble   DFDX(i,k): partial derivative of
C                                 I-th component of FCN with respect
C                                 to X(k) (output)
C       IFAIL              Int    JAC evaluation-failure indicator. 
C                                 (output)
C                                 Has always value 0 (zero) on input.
C                                 Indicates failure of JAC evaluation
C                                 and causes termination of NLSCON,
C                                 if set to a negative value on output
C
C----------------------------------------------------------------------------
C
C*   Options IOPT:
C    =============
C
C     Pos. Name   Default  Meaning
C
C       1  QSUCC  0        =0 (.FALSE.) initial call:
C                             NLSCON is not yet initialized, i.e. this
C                             is the first call for this nonlinear
C                             least squares problem.
C                             At successfull return with MODE=1,
C                             QSUCC is set to 1.
C                          =1 (.TRUE.) successive call:
C                             NLSCON is initialized already and is now
C                             called to perform one or more following
C                             Gauss-Newton-iteration steps.
C                             ATTENTION:
C                                Don't destroy the contents of
C                                IOPT(i) for 1 <= i <= 50 ,
C                                IWK(j)  for 1 <= j < NIWKFR and
C                                RWK(k)  for 1 <= k < NRWKFR.
C                                (Nevertheless, some of the options, e.g.
C                                 FCMIN, SIGMA, MPR..., can be modified
C                                 before successive calls.)
C       2  MODE   0        =0 Standard mode initial call:
C                             Return when the required accuracy for the
C                             iteration vector is reached. User defined
C                             parameters are evaluated and checked.
C                             Standard mode successive call:
C                             If NLSCON was called previously with
C                             MODE=1, it performs all remaining 
C                             iteration steps.
C                          =1 Stepwise mode:
C                             Return after one Gauss-Newton
C                             iteration step.
C       3  JACGEN 0        Method of Jacobian generation
C                          =0 Standard method is JACGEN=2
C                          =1 User supplied subroutine JAC will be 
C                             called to generate Jacobian matrix
C                          =2 Jacobian approximation by numerical
C                             differentation (see subroutine NCJAC)
C                          =3 Jacobian approximation by numerical
C                             differentation with feedback control
C                             (see subroutine NCJCF)
C       4..8               Reserved
C       9  ISCAL  0        Determines how to scale the iterate-vector:
C                          =0 The user supplied scaling vector XSCAL is
C                             used as a (componentwise) lower threshold
C                             of the current scaling vector
C                          =1 The vector XSCAL is always used as the
C                             current scaling vector
C      10                  Reserved
C      11  MPRERR 0        Print error messages
C                          =0 No output
C                          =1 Error messages
C                          =2 Warnings additionally
C                          =3 Informal messages additionally
C      12  LUERR  6        Logical unit number for error messages
C      13  MPRMON 0        Print iteration Monitor
C                          =0 No output
C                          =1 Standard output
C                          =2 Summary iteration monitor additionally
C                          =3 Detailed iteration monitor additionally
C                          =4,5,6 Outputs with increasing level addi-
C                             tional increasing information for code
C                             testing purposes. Level 6 produces
C                             in general extremely large output!
C      14  LUMON  6        Logical unit number for iteration monitor
C      15  MPRSOL 0        Print solutions
C                          =0 No output
C                          =1 Initial values and solution values
C                          =2 Intermediate iterates additionally
C      16  LUSOL  6        Logical unit number for solutions
C      17..18              Reserved
C      19  MPRTIM 0        Output level for the time monitor
C                          = 0 : no time measurement and no output
C                          = 1 : time measurement will be done and
C                                summary output will be written -
C                                regard note 4a.
C      20  LUTIM  6        Logical output unit for time monitor
C      21  QSTAT  0        Statistical Analysis of the final 
C                          least squares estimate:
C                          = 0 : Analysis will not be done
C                          = 1 : Analysis will be done,
C                                and certain results are stored
C                                to the RWK array (for details, see
C                                RWK description below)
C      22  MPRSTA 0        Printing of statistical Analysis for
C                          the final least squares estimate:
C                          = 0 : Printing will not be done
C                          = 1 : Printing will be done (and this
C                                implies QSTAT to be set to 1)
C      23..30              Reserved
C      31  NONLIN 3        Problem type specification
C                          =1 Linear problem
C                             Warning: If specified, no check will be
C                             done, if the problem is really linear, and
C                             NLSCON terminates unconditionally after
C                             one Gauss-Newton-iteration step.
C                          =2 Mildly nonlinear problem
C                          =3 Highly nonlinear problem
C                          =4 Extremely nonlinear problem
C      32  QRANK1 0        =0 (.FALSE.) Rank-1 updates by Broyden-
C                             approximation are inhibited.
C                          =1 (.TRUE.) Rank-1 updates by Broyden-
C                             approximation are allowed.
C      33..34              Reserved
C      35  QNSCAL 0        Inhibit automatic row scaling: 
C                          =0 (.FALSE.) Automatic row scaling of
C                             the linear system is activ: 
C                             Rows i=1,...,N will be divided by
C                             max j=1,...,N (abs(a(i,j))) 
C                          =1 (.TRUE.) No row scaling of the linear
C                             system. Recommended only for well row-
C                             scaled nonlinear least squares problems.
C      36  ITERM  0        Determines the iteration termination cri-
C                          terium to be chosen:
C                          =0 Iteration is terminated, if one of the
C                             stopping criteria used for ITERM=1 and
C                             ITERM=2 is satisfied (see below).
C                          =1 Iteration is terminated, if, from a
C                             statistical point of view, a resonable
C                             precision is achieved, i.e.
C                             simplified GN-correction < RTOL, and an
C                             estimate of the accuracy is available.
C                             Recommended to be used for incompatible
C                             problems.
C                          =2 Iteration is terminated, if 
C                             GN-correction < RTOL . Using this option
C                             may force 'to much' precision from
C                             the statistical point of view.
C      37                  Reserved
C      38  IBDAMP          Bounded damping strategy switch:
C                          =0 means currently always IBDAMP = off 
C                             (but may depend on the settings of other
C                              options in future versions)
C                          =1 means always IBDAMP = on 
C                          =2 means always IBDAMP = off 
C      39..45              Reserved
C      46..50              User options (see note 4b)

C     Note 3:
C         If NLSCON terminates with IERR=2 (maximum iterations)
C         or  IERR=3 (small damping factor), you may try to continue
C         the iteration by increasing NITMAX or decreasing FCMIN
C         (see RWK) and setting QSUCC to 1.
C
C     Note 4a:
C        The integrated time monitor calls the machine dependent
C        subroutine ZIBSEC to get the current time stamp in form
C        of a real number (Single precision). As delivered, this
C        subroutine always return 0.0 as time stamp value. Refer
C        to the compiler- or library manual of the FORTRAN compiler
C        which you currently use to find out how to get the current
C        time stamp on your machine.
C
C     Note 4b:
C         The user options may be interpreted by the user replacable
C         routines NCSOUT, NCFACT, NCSOLV - the distributed version
C         of NCSOUT currently uses IOPT(46) as follows:
C         0 = standard plotdata output (may be postprocessed by a user-
C             written graphical program)
C         1 = plotdata output is suitable as input to the graphical
C             package GRAZIL (based on GKS), which has been developed
C             at ZIB. 
C
C
C*   Optional INTEGER input/output in IWK:
C    =======================================
C
C     Pos. Name          Meaning
C
C      1   NITER  IN/OUT Number of Gauss-Newton-iterations
C      2                 reserved
C      3   NCORR  IN/OUT Number of corrector steps
C      4   NFCN   IN/OUT Number of FCN-evaluations
C      5   NJAC   IN/OUT Number of Jacobian generations or
C                        JAC-calls
C      6                 reserved
C      7                 reserved
C      8   NFCNJ  IN/OUT Number of FCN-evaluations for Jacobian
C                        approximation
C      9   NREJR1 IN/OUT Number of rejected Gauss-Newton iteration
C                        steps done with a rank-1 computed Jacobian
C     10..11             Reserved
C     12   IDCODE IN/OUT Output: The 8 decimal digits program identi-
C                        fication number ppppvvvv, consisting of the
C                        program code pppp and the version code vvvv.
C                        Input: If containing a negative number,
C                        it will only be overwritten by the identi-
C                        fication number, immediately followed by
C                        a return to the calling program.      
C     13..15             Reserved
C     16   NIWKFR OUT    First element of IWK which is free to be used
C                        as workspace between Gauss-Newton iterations.
C                        For standard linear solvers:  51
C     17   NRWKFR OUT    First element of RWK which is free to be used
C                        as workspace between Newton iteration steps.
C                        For standard linear solver and numerically 
C                        approximated Jacobian computed by the 
C                        expression: N*(M+4)+4*M+61 
C                        If the Jacobian is computed by a user routine
C                        JAC, subtract N in this expression.
C     18   LIWKA  OUT    Length of IWK currently required
C     19   LRWKA  OUT    Length of RWK currently required
C     20..22             Reserved
C     23   IFAIL  OUT    Set in case of failure of NCFACT (IERR=80),
C                        N2SOLV (IERR=81), FCN (IERR=82) or JAC(IERR=83)
C                        to the nonzero IFAIL value returned by the 
C                        routine indicating the failure .
C     24..30             Reserved
C     31   NITMAX IN     Maximum number of permitted iteration
C                        steps (default: 50)
C     32   IRANK  IN     Initial rank
C                        =0 : full rank N
C                        =k with 0 < k < N : deficient rank assumed
C                           for the Jacobian in the starting point
C     33   NEW    IN/OUT Count of consecutive rank-1 updates
C     34   IFCCNT INTERN Count of consecutive special steps before 
C                        convergence stop of iteration
C     35..50             Reserved
C
C*   Optional REAL input/output in RWK:
C    ====================================
C
C     Pos. Name          Meaning
C
C      1..16             Reserved
C     17   CONV   OUT    The maximum norm of the latest ordinary 
C                        respective simplified (scaled) Gauss-Newton
C                        correction.
C     18   SUMX   OUT    Natural level (not Normx of printouts)
C                        of the current iterate, i.e. Sum(DX(i)**2),
C                        where DX = scaled Gauss-Newton correction.
C     19   DLEVF  OUT    Standard level (not Normf of printouts)
C                        of the current iterate, i.e. Norm2(F(X)),
C                        where F =  nonlinear model function.
C     20   FCBND  IN     Bounded damping strategy restriction factor
C                        (Default is 10)
C     21   FCSTRT IN     Damping factor for first Gauss-Newton iteration
C                        - overrides option NONLIN, if set (see note 5)
C     22   FCMIN  IN     Minimal allowed damping factor (see note 5)
C     23   SIGMA  IN     Broyden-approximation decision parameter
C                        Required choice: SIGMA.GE.1. Increasing this
C                        parameter make it less probable that the algo-
C                        rith performs Broyden steps.
C                        Rank1 updates are inhibited, if 
C                        SIGMA.GT.1/FCMIN is set. (see note 5)
C     24                 Reserved
C     25   COND   IN     Maximum permitted subcondition for rank-
C                        decision by linear solver
C                        (Default: 1/epmach, epmach: relative
C                         machine precision) 
C     26   AJDEL  IN     Jacobian approximation without feedback:
C                        Relative pertubation for components
C                        (Default: sqrt(epmach*10), epmach: relative
C                         machine precision) 
C     27   AJMIN  IN     Jacobian approximation without feedback:
C                        Threshold value (Default: 0.0d0)
C                          The absolute pertubation for component k is
C                          computed by 
C                          DELX := AJDEL*max(abs(Xk),AJMIN)
C     28  ETADIF  IN     Jacobian approximation with feedback:
C                        Target value for relative pertubation ETA of X
C                        (Default: 1.0d-6)
C     29  ETAINI  IN     Jacobian approximation with feedback:
C                        Initial value for denominator differences
C                        (Default: 1.0d-6)
C     30                 Reserved
C     31  PREC    OUT    An estimate for the achieved relative accuracy.
C                        This number is only available, if IERR=0 or 
C                        IERR=1 and an estimate for the incompatibility
C                        factor kappa (SKAP=RWK(32)) can be made. If no
C                        meaningful information is at hand, a -1.0 is
C                        stored.
C     32  SKAP    OUT    An estimate of the incompatibility factor of
C                        the given nonlinear least squares problem.
C                        This number is only available, if IERR=0 or 
C                        IERR=1 and a certain further condition holds.
C                        If no meaningful information is at hand, a -1.0
C                        is stored.
C                        A small value of SKAP indicates that the given
C                        measurement data can be well fitted by the
C                        model, whereas a value SKAP.GE.0.5 may give a 
C                        hint to an insufficient modelling of the
C                        experiment from which the measurements
C                        originate.
C     33..49             Reserved
C     50  SIGMA2  OUT    Holds the estimated variance of the residual
C                        on final exit of NLSCON, if IOPT(21)=1 is set.
C     51..50+N
C         XL(N)          Holds the left bounds of the final parameters
C                        confidence intervals; 
C                        on final exit of NLSCON, if IOPT(21)=1 is set.
C     50+N+1..50+2*N
C         XR(N)          Holds the right bounds of the final parameters
C                        confidence intervals; 
C                        on final exit of NLSCON, if IOPT(21)=1 is set.
C     50+2*N+1..50+2*N+N*N         
C        VCV(N,N) OUT    Holds the columnwise stored correlation-matrix
C                        on final exit of NLSCON, if IOPT(21)=1 is set.
C
C     Note 5:
C       The default values of the internal parameters may be obtained
C       from the monitor output with at least IOPT field MPRMON set to 2
C       and by initializing the corresponding RWK-fields to zero. 
C
C*   Error messages:
C    ===============
C
C      1    Termination at stationary point (rank deficient Jacobian
C           and termination criterion fulfilled)
C      2    Termination after NITMAX iterations ( as indicated by
C           input parameter NITMAX=IWK(31) )
C      3    Termination, since damping factor became to small and
C           Jacobian rank was already reduced to zero
C     10    Integer or real workspace too small
C     20    Bad or inconsistent input to one or more of the 
C           dimensional parameters N,M and MFIT
C     21    Nonpositive value for RTOL supplied
C     22    Negative scaling value via vector XSCAL supplied
C     30    One or more fields specified in IOPT are invalid
C           (for more information, see error-printout)
C     80    Error signalled by linear solver routine NCFACT,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C           (not used with standard routine NCFACT)
C     81    Error signalled by linear solver routine NCFIT,
C           for more detailed information see IFAIL-value
C           stored to IWK(23)
C           (not used with standard routine NCFIT)
C     82    Error signalled by user routine FCN (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C     83    Error signalled by user routine JAC (Nonzero value
C           returned via IFAIL-flag; stored to IWK(23) )
C     180,182,183
C           see error codes 80,82,83, but the failure of NCFACT, FCN
C           or JAC occured when preparing the call of the statistics
C           subroutine STACON for the final iterate of NLSCON.
C
C     Note 6 : in case of failure:
C        -    use non-standard options
C        -    use another initial guess
C        -    or reformulate model
C        -    or turn to general optimization routine
C
*/

typedef void (*FcnNlscon)(int* n, int* m, int* mcon,
                           double* x, double* f, int* ifail);

typedef void (*JacNlscon)(int* n, int* m, int* mcon,
                           double* x, double* dfdx, int* ifail);


typedef struct 
  {
    int       n;
    int       m;
    int       mFit;
    FcnNlscon fcn;
    JacNlscon jac;
    double    *x, *xScal;  /* n unknowns */
    double    *fi, *fScal; /* mFit data  */
    double    rTol;        /* rTol w.r.t. unknown x */
    int       iOpt[50];
    int       iErr;
    int       lIwk;        /* lIwk = n + 52 */
    int       *iwk;
    int       lRwk;        /* lRwk = 2*(m+n)*n + 8*m + 10*n + max(m,n) + 104 */
    double    *rwk;
  } 
nlscon_state_t;


/* NLSCON API */

extern void nlscon_(
            int*      n,
            int*      m,
            int*      mFit,
            FcnNlscon fcn,
            JacNlscon jac,
            double*   x,
            double*   xScal,
            double*   fi,
            double*   fScal,
            double*   rTol,
            int*      iOpt,
            int*      iErr,
            int*      lIwk,
            int*      iwk,
            int*      lRwk,
            double*   rwk);

