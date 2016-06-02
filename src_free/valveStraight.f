      subroutine valveStraight(x, prmt)
!     Calculation of valve opening is seperated from outp (Compared
!     with original version.)
!     x and prmt is necessary to be included because:
!     x is used to calculate the current timing;
!     prmt(3) is the current step width.
!     CV() is already a global variable declared in data.f.

      use data
      implicit none

      double precision x,prmt
      dimension prmt(8)

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

!     Straight line valve opening
      IF (XT.GE.XX(1) .AND. XT.LT.XX(2)) THEN
         XK1=(XT-XX(1))/(XX(2)-XX(1))
         IF(XK1.LT.0.5) THEN
            B=2.0D0*PI*XK1
         ELSE
            B=2.0D0*PI-2.0D0*PI*XK1
         END IF
         CV(1)=CV10*B
      END IF
      IF (XT.GE.XX(2)) THEN
         CV(1)=0.0D0
      ENDIF

!     Straight line valve opening
      IF (XT.GE.XX(3) .AND. XT.LT.XX(4)) THEN
         XK2=(XT-XX(3))/(XX(4)-XX(3))
         IF(XK2.LT.0.5) THEN
            B=2.0D0*PI*XK2
         ELSE
            B=2.0D0*PI-2.0D0*PI*XK2
         END IF
         CV(2)=CV20*B
      END IF
      IF (XT.GE.XX(4)) THEN
         CV(2)=0.0D0
      ENDIF
      end subroutine
