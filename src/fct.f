	SUBROUTINE FCT(X,Y,DERY)
!     This subroutine is a "MUST" for Runge-Kutta subroutine.
!     Look http://media.ibm1130.org/1130-006-ocr.pdf for detailed description.
!     Here, the right hand side of the differential equation, namely y'(x)
!     is calculated.

	use data
	implicit none
	double precision y, dery, prmt, x
	DIMENSION Y(4), dery(4), PRMT(8)

!       The original version of Prof. Matsubara does not calculate
!       valve opening here. It may cause slight difference in result.
	call valve(x, prmt)

	V2=0.5*VT*(1.0-COS(OMEGA*X-PHAI))
	V1=VT-V2
	P1=(C1+Y(1)-Y(2)-Y(3))*R/(VRH/TRH+V1/TH)
	P2=(C2+Y(3))*R/(VRE/TRE+V2/TE)
	PP1=PH*PH-P1*P1
	PP2=P1*P1-PL*PL
!       PPX1=SIGN(1.0,PP1)*SQRT(ABS(PP1))
!       PPX2=SIGN(1.0,PP2)*SQRT(ABS(PP2))
	PPX1=SQRT(ABS(PP1))
	PPX2=SQRT(ABS(PP2))
	DERY(1)=CV(1)*PPX1
	DERY(2)=CV(2)*PPX2
!       DERY(3)=CV(3)*(P1*P1-P2*P2)
	DERY(3)=CV(3)*(P1-P2)

!       Dummy space holder for future displaser calculation.
	DERY(4)=0.0D0

        RETURN
        END SUBROUTINE
