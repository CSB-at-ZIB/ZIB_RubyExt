      SUBROUTINE ZIBCONST(EPMACH,SMALL)
      DOUBLE PRECISION EPMACH,SMALL
C
C*********************************************************************
C  Set Approximations to machine constants.
C  This routine is machine dependent.
C*********************************************************************
C
C  Output parameters
C    EPMACH    DOUBLE     Relative machine precision
C    SMALL     DOUBLE     SQRT of smallest positive machine number
C 
C  The proposed values should work on Intel compatible cpus,
C  PowerPCs and Sun sparcs
C
C*********************************************************************
C
      EPMACH = 2.22D-16
      SMALL  = 1.0D-145
C     EPMACH = 1.0D-17
C     SMALL  = 1.0D-150
      RETURN
      END
