	SUBROUTINE FCT(X,Y,DERY)
!     This subroutine is a "MUST" for Runge-Kutta subroutine.
!     Look http://media.ibm1130.org/1130-006-ocr.pdf for detailed description.
!     Here, the right hand side of the differential equation, namely y'(x)
!     is calculated.

	use data
	implicit none
	double precision y, dery, prmt, x
	DIMENSION Y(5), dery(5), PRMT(8)

	double precision pressureAssist
	double precision area

!       The original version of Prof. Matsubara does not calculate
!       valve opening here. It may cause slight difference in result.

!       valve() is the original version (Triangle opening)
!       valveTrapezoid is the new trapezoid shape opening
!	call valve(x, prmt)
	call valveTrapezoid(x, prmt)

!       Pressure in assist room is calculated by a function.
	pa = pressureAssist(x,f)

!       Original version of volume calculation:
!	V2=0.5*VT*(1.0-COS(OMEGA*X-PHAI))
!	V1=VT-V2

!       V1 - Compression room volume. (cm^3)
!       V2 - Expansion room volume. (cm^3)
!       A gas driven calculation: ( y(4)*100 means convert m -> cm )
!       VT/4.0 because the current stroke is 40 mm, thus 4 cm.
	if (drive_mode == '1') then
	   V2=0.5*VT*(1.0-COS(OMEGA*X-PHAI))
	   V1=VT-V2
	else if (drive_mode == '2') then
	   Area = VT/4.0d0
!	One can add dead volume to V2, which will offset the volume
!	curve. But the final result won't change much.
	   V2=0.5d0*VT - y(4)*100*Area
	   V1=VT-V2
	end if
	P1=(C1+Y(1)-Y(2)-Y(3))*R/(VRH/TRH+V1/TH)
	P2=(C2+Y(3))*R/(VRE/TRE+V2/TE)
	PP1=PH*PH-P1*P1
	PP2=P1*P1-PL*PL

!	It's important to allow reverse flow across the valve.
!	Otherwise pressure overshooting (ex. Pressure over the PH)
!	might occur.
        PPX1=SIGN(1.0D0,PP1)*SQRT(ABS(PP1))
        PPX2=SIGN(1.0D0,PP2)*SQRT(ABS(PP2))

!       DERY(1) -  dm1 / dt - mass flow through rotary valve (High Pressure).
!       DERY(2) -  dm2 / dt - mass flow through rotary valve (Low Pressure).
	DERY(1)=CV(1)*PPX1
	DERY(2)=CV(2)*PPX2

!       DERY(3) - dm / dt - mass flow through regenerator
!       DERY(3)=CV(3)*(P1*P1-P2*P2)
	DERY(3)=CV(3)*(P1-P2)

!       Dummy space holder for future displaser calculation.
!       The displacement velocity derived from summed force.
!       y(4) - position
!       y(5) = dery(4) - velocity
	dery(4)=y(5)

!       Pressure - MPa
!       Area - m^2
!       Mass - kg
!       Displacement - m

!       displacerMass - kg
!       dampCoef - Resist coefficient on velocity term.
!       springCoef - Resist coefficient on position term.
!       y'' = deltaP*A/m - veloResit/m * y' + g
!       
!       Scotch-yoke cross area: 0.25 * 3.14 * 0.045^2 = 1.589e-3
!	However, note that in previous exp, the 45mm was used
!	to drive 2 displacers, instead of 1.
!	So the velocity result here will be significantly higher
!	since only one displacer.

	dery(5)=((pa-p1)*1000000*0.00159D0) / displacerMass -
     +  (dampCoef / displacerMass) * y(5) - springCoef * y(4) +
     +  9.80D0

!	If the displacer hit the limitation...
	if (((y(4).ge.0.02D0).or.(y(4).le.-0.02D0))) then
!       If the velocity direction and displacement is the same,
!       Then it should be stopped.
	   if (y(4)*dery(4).gt.0.0D0) then
	   dery(4) = 0.0D0
	   end if
	end if

        RETURN
        END SUBROUTINE
