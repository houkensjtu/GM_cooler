	PROGRAM GMA
!       A simulation program for simplified GM cycle cryocooler
!       originally created by Y. Matsubara.
!       Look http://ci.nii.ac.jp/naid/130004457265 ( the original paper)
!       for further informtion.

!       Variable declaration is now seperated to a external module.
!       In order to avoid conflicts of data type and ambiguous 
!       definition.
	use data
	implicit none

!       These variables are used by the Runge-Kutta subroutine.
	double precision y, dery, prmt, aux, x
	DIMENSION Y(4),DERY(4),AUX(8,4),PRMT(8)
	integer ndim, ihlf

	OPEN(8,FILE='GMout.txt')
	OPEN(9,FILE='GM.txt')
	OPEN(3,FILE='IN.TXT')
	READ(3,*)VRD,VT,PH,PL,F,FIE
	READ(3,*)B1,B2,B3,MODE
	READ(3,*)X1,X2,X3,X4
	READ(3,*)PRMT(2),PRMT(3)
	READ(3,*)R,TH,TE
!       PHAI0 : Delay angle
!       JD    : Data pitch
	READ(3,*)PHAI0,JD

!       Set global variables as parameters.
	call parSetting(y, dery, prmt, ndim)

	CALL RKGS(PRMT,Y,DERY,NDIM,IHLF,AUX)

!       Writing out data to GM.txt
	JJ = log_step
!	Old F77 style line continuation.
	WRITE(9,101)'Deg.','Ve','P1','P2','dmr/dt','dm1/dt','dm2/dt','Cv1'
!	Note that the second line should start from behind the "WRITE".
     +              ,'Cv2','Pos','P1-P2'
101	FORMAT(11A12)
!	JD=100
	DO J=1,JJ,JD
	K=360.0*J/JJ
	DP=D(2,J)-D(3,J)
	WRITE(9,100)K,(D(I,J),I=1,9),DP
	END DO
100	FORMAT(I12,10F12.3)

!       Below is the output part.
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
	
	! Hot end.
	WE1=0.5*(D(2,1)+D(2,JJ))*(D(1,1)-D(1,JJ))
	! Cold end.
	WE2=0.5*(D(3,1)+D(3,JJ))*(D(1,1)-D(1,JJ))

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

	! Work flow at Ph
	WF1=WF1+D(5,J)*DLOG(PH/D(2,J))
	WU1=WU1+D(5,J)/D(2,J)*DLOG(PH/D(2,J))
	PM=PM+D(2,J)
	WF2=WF2+D(6,J)*DLOG(D(2,J)/PL)

	! Pl
	WU2=WU2+D(6,J)/D(2,J)*DLOG(D(2,J)/PL)

	! REG
	WF0=WF0+D(4,J)*DLOG(D(2,J)/D(3,J))

	! Work input
	WF10=WF10+D(5,J)*DLOG(PH/PL)
	WF20=WF20+D(6,J)*DLOG(PH/PL)
3	CONTINUE
	AMR=AMR/JJ
	PM=PM/JJ

!	FM1=FM1/JJ*2.0 <- must be removed 2.0

	FM1=FM1/JJ			
	FM2=FM2/JJ
	HR=AMR*CP*(TH-TE)*FIE
!	HR=F*AMEMAX*CP*(TH-TE)*FIE

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

	!Actual input power using comp efficiency of 55%
	WACT=WF10/0.55

	!CHANGED 2006 FROM AM0 TO AMR
	WRITE(*,200)' m,We(J/s),,Q=',AMR,WE1,WE2,Q

200	FORMAT(A,4F10.3)
C--------------------------------------
	WRITE(8,201)' Wa,We,Qe=',WE1,WE2,Q
	WRITE(8,202)' WB1,WB2,Wr=',WF1,WF2,WF0
	WRITE(*,202)' WB1,WB2,Wr,Wu1=',WF1,WF2,WF0,WU1
	WRITE(8,203)' <W1>,<W2>,Wact=',WF10,WF20,WACT
201	FORMAT(A,4F14.5)
202	FORMAT(A,4F14.5)
203	FORMAT(A,3F12.5)

!       No assumption.
	CARACT=Q/WACT/(TE/(TH-TE))*100.0 

!       Assume ideal compressor
	CARNOT=Q/WF10/(TE/(TH-TE))*100.0

!       Assume ideal regenerator
	CMAX=WE2/WF10/(TE/(TH-TE))*100.0

!       Assume ideal reg and valve
	CMAX2=WE2/(WF10-WF1-WF2)/(TE/(TH-TE))*100.0
	WRITE(8,109)' dm1/dt,dm2/dt=',FM1,FM2
	WRITE(*,109)' FM1,2=',FM1,FM2
	WRITE(8,108)' %Carnot(1,2,3,4)=',CARACT,CARNOT,CMAX,CMAX2
	WRITE(*,108)' %Carnot=',CARACT,CARNOT,CMAX,CMAX2
108	FORMAT(A,4F10.2)
109	FORMAT(A,2F10.4) 

	END PROGRAM



