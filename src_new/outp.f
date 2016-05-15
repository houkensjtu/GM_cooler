      SUBROUTINE OUTP(X,Y,DERY,IHLF,NDIM,PRMT)
!     This subroutine is a "MUST" for Runge-Kutta subroutine.
!     Look http://media.ibm1130.org/1130-006-ocr.pdf for detailed description.
!     However, instead of writing out data, this subroutine is mainly
!     aimed at updating the value of Valve opening - CV(3).

      use data
      implicit none
      double precision y, dery, prmt, aux, x
      DIMENSION Y(3),DERY(3),AUX(8,3),PRMT(8)
      integer ndim, ihlf

      IF(X.EQ.0.) THEN
         PI=3.141592
         log_step = 0
         N=1
         FF=1./F
         XT=0.
         PRMT(8)=0.
         CV(1)=0.
         CV(2)=0.
      ENDIF

      NX=X*1000
      NF=FF*1000
      IF(NX.EQ.NF) THEN
!     WRITE(*,*)N,X,VE0,VI0
         N=N+1
         FF=N/F
         XT=0.
         XK1=0.
         XK2=0.
         XK3=0.
         XK4=0.
      ENDIF
      XT=XT+PRMT(3)

      IF (XT.GE.XX(1) .AND. XT.LT.XX(2)) THEN
         XK1=(XT-XX(1))/(XX(2)-XX(1))
!     Following part was changed to square opening
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
!--------------------------------------
C     AH=ABS(1.0-2.0*XK1)
C     IF(AH.GT.1.0) AH=1.0
C     S=2.0*ACOS(AH)
C     A=S-SIN(S)
C     CV(1)=CV10*A
C     ENDIF
!--------------------------------------

      IF (XT.GE.XX(2)) THEN
         CV(1)=0.
      ENDIF
      IF (XT.GE.XX(3) .AND. XT.LT.XX(4)) THEN
         XK2=(XT-XX(3))/(XX(4)-XX(3))
!--------------------------------------
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

!--------------------------------------
C     AH=ABS(1.0-2.0*XK2)
C     IF(AH.GT.1.0) AH=1.0
C     S=2.0*ACOS(AH)
C     A=S-SIN(S)
C     CV(2)=CV20*A
C     ENDIF
!--------------------------------------

      IF (XT.GE.XX(4)) THEN
         CV(2)=0.
      ENDIF
!-------------------------------------
!     Prmt(7) = Xend - 0.5*(Full cycle period)
!     Means to output only the final cycle to GM.txt
      IF(X.GE.PRMT(7)) THEN
         log_step = log_step + 1
         D(1,log_step)=V2

!        Pressure
         D(2,log_step)=P1

!        Pressure
         D(3,log_step)=P2

!        mass flow rate at regenerator
         D(4,log_step)=DERY(3)

!        mass flow rate at Valve1
         D(5,log_step)=DERY(1)

!         mass flow rate at Valve2
         D(6,log_step)=DERY(2)

!        Valve 1 opening
         D(7,log_step)=CV(1)/21.2

!        Valve 2 opening
         D(8,log_step)=CV(2)/21.2
      ENDIF
!-------------------------------------
      RETURN
      END SUBROUTINE
