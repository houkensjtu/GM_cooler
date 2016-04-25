	PROGRAM GMA
C-----------------------------------------
c
	use data
	implicit none
	double precision y, dery, prmt, aux, x
	DIMENSION Y(3),DERY(3),AUX(8,3),PRMT(8)
	integer ndim, ihlf

	OPEN(8,FILE='GMout.txt')
	OPEN(9,FILE='GM.txt')
	OPEN(3,FILE='IN.TXT')
	REWIND 3
	READ(3,*)VRD,VT,PH,PL,F,FIE
	READ(3,*)B1,B2,B3,MODE				 !Vh, Vl, reg, 0 for RV
	CV(1)=B1*21.2
	CV(2)=B2*21.2
	CV(3)=B3
	READ(3,*)X1,X2,X3,X4                     !Degree(<360)
	READ(3,*)PRMT(2),PRMT(3)
	PRMT(1)=0.0
	PRMT(4)=0.0001
	READ(3,*)R,TH,TE
	READ(3,*)PHAI0,JD		!DELAY ANGLE(DEG),DATA PICTH FOR FIG.
	PI=3.141592
	OMEGA=2.0*PI*F
	PHAI=PHAI0*PI/180.0
	CP=5.2
	TR=0.5*(TH+TE)
	P0=0.5*(PH+PL)
        TRH=(3.0*TH+TE)/4.0
	TRE=(3.0*TE+TH)/4.0
	VRH=0.5*VRD
	VRE=0.5*VRD
	PM=(PH+PL)/2.0
	V20=0.5*VT*(1.0-COS(-PHAI))
	V10=VT-V20
	C1=PL/R*(VRH/TRH+V10/TH)
	C2=PL/R*(VRE/TRE+V20/TE)
	CV10=CV(1)
	CV20=CV(2)
	XX(1)=X1/F/360.
	XX(2)=X2/F/360.
	XX(3)=X3/F/360.
	XX(4)=X4/F/360.
	PRMT(7)=PRMT(2)-1./F
C--- R=2.077,(for He [MPa,cm3/g/K] ---
	Y(1)=0.0
	Y(2)=0.0
	Y(3)=0.0
	DERY(1)=0.3
	DERY(2)=0.3
	DERY(3)=0.4
	NDIM=3
C	 
C      CALL RKGS(PRMT,Y,DERY,NDIM,IHLF,AOX)
       CALL rkgs(PRMT,Y,DERY,NDIM,IHLF,aux)
	JJ=J
	WRITE(9,101)'Deg.',' Ve',' P1',' P2',' dmr/dt',' dm1/dt',' dm2/dt'
     + ,' Cv1',' Cv2',' P1-P2'
101	FORMAT(10A12)
C==================================
	AMR=0.
	AME=0.
	FM1=0.0
	FM2=0.0
	WF0=0.
        WF1=0.
	WU1=0.
	PM=0.
      WF2=0.
      WF10=0.
      WF20=0.
      WE1=0.5*(D(2,1)+D(2,JJ))*(D(1,1)-D(1,JJ))		!HOT END
	WE2=0.5*(D(3,1)+D(3,JJ))*(D(1,1)-D(1,JJ))		!COLD END
	DO 3 J=1,JJ
	FM1=FM1+D(5,J)
	FM2=FM2+D(6,J)
	AME=D(1,J)*D(3,J)/R/TE
	IF(AME.GT.AMEMAX) AMEMAX=AME
	AMR=AMR+ABS(D(4,J))
	IF(D(4,J).GT.AM0) AM0=D(4,J)
	IF(J.GT.1) THEN
	WE1=WE1+0.5*(D(2,J)+D(2,J-1))*(D(1,J)-D(1,J-1))
	WE2=WE2+0.5*(D(3,J)+D(3,J-1))*(D(1,J)-D(1,J-1))
	ENDIF
	WF1=WF1+D(5,J)*DLOG(PH/D(2,J))         !Work flow at Ph
	WU1=WU1+D(5,J)/D(2,J)*DLOG(PH/D(2,J))
	PM=PM+D(2,J)
	WF2=WF2+D(6,J)*DLOG(D(2,J)/PL)
	WU2=WU2+D(6,J)/D(2,J)*DLOG(D(2,J)/PL)         !             Pl
	WF0=WF0+D(4,J)*DLOG(D(2,J)/D(3,J))     !             REG
	WF10=WF10+D(5,J)*DLOG(PH/PL)           !Work input
	WF20=WF20+D(6,J)*DLOG(PH/PL)
3	CONTINUE
	 AMR=AMR/JJ
	PM=PM/JJ
c	FM1=FM1/JJ*2.0			!must be removed 2.0
	FM1=FM1/JJ			
	FM2=FM2/JJ
	HR=AMR*CP*(TH-TE)*FIE
C	HR=F*AMEMAX*CP*(TH-TE)*FIE
	WEE=WE2/AM0
      WE1=WE1*F
	WE2=WE2*F
	Q=WE2-HR
	WF1=WF1*R*TH/JJ
	WU1=WU1*R*TH*PM/JJ
	WF2=WF2*R*TH/JJ
	WF0=WF0*R*TR/JJ
	WF10=WF10*R*TH/JJ
	WF20=WF20*R*TH/JJ
	WACT=WF10/0.55			!Actual input power using comp efficiency of 55%
	WRITE(*,200)' m,We(J/s),,Q=',AMR,WE1,WE2,Q	!CHANGED 2006 FROM AM0 TO AMR
200	FORMAT(A,4F10.3)
C--------------------------------------
	WRITE(8,201)' Wa,We,Qe=',WE1,WE2,Q
      WRITE(8,202)' WB1,WB2,Wr=',WF1,WF2,WF0
      WRITE(*,202)' WB1,WB2,Wr,Wu1=',WF1,WF2,WF0,WU1
      WRITE(8,203)' <W1>,<W2>,Wact=',WF10,WF20,WACT
201   FORMAT(A,4F14.5)
202   FORMAT(A,4F14.5)
203   FORMAT(A,3F12.5)
	CARACT=Q/WACT/(TE/(TH-TE))*100.0			!No asumptions
      CARNOT=Q/WF10/(TE/(TH-TE))*100.0			!Assume ideal compressor
	CMAX=WE2/WF10/(TE/(TH-TE))*100.0			!Assume ideal regenerator
	CMAX2=WE2/(WF10-WF1-WF2)/(TE/(TH-TE))*100.0	!Assume ideal reg and valve
	WRITE(8,109)' dm1/dt,dm2/dt=',FM1,FM2
	WRITE(*,109)' FM1,2=',FM1,FM2
	WRITE(8,108)' %Carnot(1,2,3,4)=',CARACT,CARNOT,CMAX,CMAX2
	WRITE(*,108)' %Carnot=',CARACT,CARNOT,CMAX,CMAX2
108	FORMAT(A,4F10.2)
109   FORMAT(A,2F10.4) 
C===========================
C	JD=100
	DO J=1,JJ,JD
	K=360.0*J/JJ
	DP=D(2,J)-D(3,J)
	WRITE(9,100)K,(D(I,J),I=1,8),DP
	END DO
100	FORMAT(I12,9F12.3)
	END
C----------------------------------
      SUBROUTINE FCT(X,Y,DERY)
	use data
	implicit none
	double precision x, y, dery
	dimension y(3), dery(3)

	V2=0.5*VT*(1.0-COS(OMEGA*X-PHAI))
	V1=VT-V2
	P1=(C1+Y(1)-Y(2)-Y(3))*R/(VRH/TRH+V1/TH)
	P2=(C2+Y(3))*R/(VRE/TRE+V2/TE)
	PP1=PH*PH-P1*P1
	PP2=P1*P1-PL*PL
C	PPX1=SIGN(1.0,PP1)*SQRT(ABS(PP1))
C	PPX2=SIGN(1.0,PP2)*SQRT(ABS(PP2))
	PPX1=SQRT(ABS(PP1))
	PPX2=SQRT(ABS(PP2))
	DERY(1)=CV(1)*PPX1
	DERY(2)=CV(2)*PPX2
c	DERY(3)=CV(3)*(P1*P1-P2*P2)
	DERY(3)=CV(3)*(P1-P2)
      RETURN
      END
C-----------------------------------
      SUBROUTINE OUTP(X,Y,DERY,IHLF,NDIM,PRMT)
	use data
	implicit none
	double precision x, y, dery, prmt
	dimension y(3), dery(3), prmt(8)
	integer ihlf, ndim
	IF(X.EQ.0.) THEN
	PI=3.141592
	J=0
	N=1
	FF=1./F
	XT=0.
	PRMT(8)=0.
	CV(1)=0.
	CV(2)=0.
	ENDIF
C----------------------------------------------------------------
	NX=X*1000
	NF=FF*1000
	IF(NX.EQ.NF) THEN
C	WRITE(*,*)N,X,VE0,VI0
	N=N+1
	FF=N/F
	XT=0.
	XK1=0.
	XK2=0.
	XK3=0.
	XK4=0.
	ENDIF
	  XT=XT+PRMT(3)
C---------------
	IF (XT.GE.XX(1) .AND. XT.LT.XX(2)) THEN
	XK1=(XT-XX(1))/(XX(2)-XX(1))
c--------------------------------------<<Following part was changed to square opening
	IF(XK1.LT.0.5) THEN
	C=ABS(1.0-4.0*XK1)
	S=2.0*ACOS(C)
	B=PI-0.5*(S-SIN(S))
	IF(XK1.LT.0.25) B=0.5*(S-SIN(S))
	ELSE
	C=ABS(3.0-4.0*XK1)
	S=2.0*ACOS(C)
	B=0.5*(S-SIN(S))
	IF(XK1.LT.0.75) B=PI-0.5*(S-SIN(S))
	END IF
	IF(MODE.EQ.1) B=1.0
	CV(1)=CV10*B
	END IF
C--------------------------------------
C	AH=ABS(1.0-2.0*XK1)
C	IF(AH.GT.1.0) AH=1.0
C	S=2.0*ACOS(AH)
C	A=S-SIN(S)
C	CV(1)=CV10*A
C	ENDIF
	IF (XT.GE.XX(2)) THEN
	CV(1)=0.
	ENDIF
	IF (XT.GE.XX(3) .AND. XT.LT.XX(4)) THEN
	XK2=(XT-XX(3))/(XX(4)-XX(3))
c--------------------------------------
	IF(XK2.LT.0.5) THEN
	C=ABS(1.0-4.0*XK2)
	S=2.0*ACOS(C)
	B=PI-0.5*(S-SIN(S))
	IF(XK2.LT.0.25) B=0.5*(S-SIN(S))
	ELSE
	C=ABS(3.0-4.0*XK2)
	S=2.0*ACOS(C)
	B=0.5*(S-SIN(S))
	IF(XK2.LT.0.75) B=PI-0.5*(S-SIN(S))
	END IF
	IF(MODE.EQ.1) B=1.0
	CV(2)=CV20*B
	END IF
C--------------------------------------
C	AH=ABS(1.0-2.0*XK2)
C	IF(AH.GT.1.0) AH=1.0
C	S=2.0*ACOS(AH)
C	A=S-SIN(S)
C	CV(2)=CV20*A
C      ENDIF
	IF (XT.GE.XX(4)) THEN
	CV(2)=0.
	ENDIF
C---------------------
C	IF(X.GE.PRMT(7)) THEN
	IF(X.ge.0.0) THEN
C------------------------------------
        J=J+1
	D(1,J)=V2
	D(2,J)=P1           !Press 
	D(3,J)=P2           !Press 
	D(4,J)=DERY(3)         !mass flow rate at regenerator
	D(5,J)=DERY(1)         !mass flow rate at V1
	D(6,J)=DERY(2)         !mass flow rate at V2
	D(7,J)=CV(1)/21.2           !VALVE 
	D(8,J)=CV(2)/21.2           !VALVE 
	ENDIF
C---------------------
	RETURN
	END

	
      SUBROUTINE RKGS(PRMT,Y,DERY,NDIM,IHLF,AUX)
C
C
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

