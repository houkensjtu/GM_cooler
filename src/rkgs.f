      SUBROUTINE RKGS(PRMT,Y,DERY,NDIM,IHLF,AUX)
!     A 4th order Runge-Kutta solver for DE (Differential Equation) matrix.
!     All the calculation program is unchanged from original version.
!     Look http://media.ibm1130.org/1130-006-ocr.pdf for detailed description.
!     Variable declaration is changed to "implicit none" and all variables
!     are defined before using.
      implicit none
      double precision a, b, c, xend, h, aj, bj, cj, r1, r2, delt
      integer irec, istep, itest, imod, iend, i, j, ndim, ihlf
      double precision y, dery, aux, prmt, x
      DIMENSION Y(3),DERY(3),AUX(8,3),A(4),B(4),C(4),PRMT(8)

      DO 1 I=1,NDIM
    1 AUX(8,I)=.06666667*DERY(I)
      X=PRMT(1)
      XEND=PRMT(2)
      H=PRMT(3)
      PRMT(5)=0.
      CALL FCT(X,Y,DERY)
C
C     ERROR TEST
      IF(H*(XEND-X))38,37,2
C
C     PREPARATIONS FOR RUNGE-KUTTA METHOD
    2 A(1)=.5
      A(2)=.2928932
      A(3)=1.707107
      A(4)=.1666667
      B(1)=2.
      B(2)=1.
      B(3)=1.
      B(4)=2.
      C(1)=.5
      C(2)=.2928932
      C(3)=1.707107
      C(4)=.5
C
C     PREPARATIONS OF FIRST RUNGE-KUTTA STEP
      DO 3 I=1,NDIM
      AUX(1,I)=Y(I)
      AUX(2,I)=DERY(I)
      AUX(3,I)=0.
    3 AUX(6,I)=0.
      IREC=0
      H=H+H
      IHLF=-1
      ISTEP=0
      IEND=0
C
C
C     START OF A RUNGE-KUTTA STEP
    4 IF((X+H-XEND)*H)7,6,5
    5 H=XEND-X
    6 IEND=1
C
C     RECORDING OF INITIAL VALUES OF THIS STEP
    7 CALL OUTP(X,Y,DERY,IREC,NDIM,PRMT)
      IF(PRMT(5))40,8,40
    8 ITEST=0
    9 ISTEP=ISTEP+1
C
C
C     START OF INNERMOST RUNGE-KUTTA LOOP
      J=1
   10 AJ=A(J)
      BJ=B(J)
      CJ=C(J)
      DO 11 I=1,NDIM
      R1=H*DERY(I)
      R2=AJ*(R1-BJ*AUX(6,I))
      Y(I)=Y(I)+R2
      R2=R2+R2+R2
   11 AUX(6,I)=AUX(6,I)+R2-CJ*R1
      IF(J-4)12,15,15
   12 J=J+1
      IF(J-3)13,14,13
   13 X=X+.5*H
   14 CALL FCT(X,Y,DERY)
      GOTO 10
C     END OF INNERMOST RUNGE-KUTTA LOOP
C
C
C     TEST OF ACCURACY
   15 IF(ITEST)16,16,20
C
C     IN CASE ITEST=0 THERE IS NO POSSIBILITY FOR TESTING OF ACCURACY
   16 DO 17 I=1,NDIM
   17 AUX(4,I)=Y(I)
      ITEST=1
      ISTEP=ISTEP+ISTEP-2
   18 IHLF=IHLF+1
      X=X-H
      H=.5*H
      DO 19 I=1,NDIM
      Y(I)=AUX(1,I)
      DERY(I)=AUX(2,I)
   19 AUX(6,I)=AUX(3,I)
      GOTO 9
C
C     IN CASE ITEST=1 TESTING OF ACCURACY IS POSSIBLE
   20 IMOD=ISTEP/2
      IF(ISTEP-IMOD-IMOD)21,23,21
   21 CALL FCT(X,Y,DERY)
      DO 22 I=1,NDIM
      AUX(5,I)=Y(I)
   22 AUX(7,I)=DERY(I)
      GOTO 9
C
C     COMPUTATION OF TEST VALUE DELT
   23 DELT=0.
      DO 24 I=1,NDIM
   24 DELT=DELT+AUX(8,I)*ABS(AUX(4,I)-Y(I))
      IF(DELT-PRMT(4))28,28,25
C
C     ERROR IS TOO GREAT
   25 IF(IHLF-10)26,36,36
   26 DO 27 I=1,NDIM
   27 AUX(4,I)=AUX(5,I)
      ISTEP=ISTEP+ISTEP-4
      X=X-H
      IEND=0
      GOTO 18
C
C     RESULT VALUES ARE GOOD
   28 CALL FCT(X,Y,DERY)
      DO 29 I=1,NDIM
      AUX(1,I)=Y(I)
      AUX(2,I)=DERY(I)
      AUX(3,I)=AUX(6,I)
      Y(I)=AUX(5,I)
   29 DERY(I)=AUX(7,I)
      CALL OUTP(X-H,Y,DERY,IHLF,NDIM,PRMT)
      IF(PRMT(5))40,30,40
   30 DO 31 I=1,NDIM
      Y(I)=AUX(1,I)
   31 DERY(I)=AUX(2,I)
      IREC=IHLF
      IF(IEND)32,32,39
C
C     INCREMENT GETS DOUBLED
   32 IHLF=IHLF-1
      ISTEP=ISTEP/2
      H=H+H
      IF(IHLF)4,33,33
   33 IMOD=ISTEP/2
      IF(ISTEP-IMOD-IMOD)4,34,4
   34 IF(DELT-.02*PRMT(4))35,35,4
   35 IHLF=IHLF-1
      ISTEP=ISTEP/2
      H=H+H
      GOTO 4
C
C
C     RETURNS TO CALLING PROGRAM
   36 IHLF=11
      CALL FCT(X,Y,DERY)
      GOTO 39
   37 IHLF=12
      GOTO 39
   38 IHLF=13
   39 CALL OUTP(X,Y,DERY,IHLF,NDIM,PRMT)
   40 RETURN
      END
C
C=======================================================================
C

