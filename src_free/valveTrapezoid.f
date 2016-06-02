      subroutine valveTrapezoid(x, prmt)
!     Calculation of valve opening is seperated from outp (Compared
!     with original version.)
!     x and prmt is necessary to be included because:
!     x is used to calculate the current timing;
!     prmt(3) is the current step width.
!     CV() is already a global variable declared in data.f.
!     
!     The original version is in valve.f, which only allows triangle
!     opening. In this new subroutine, 4 new angle parameters
!     x5~x8 were added and the opening shape is changed to trapezoid.

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
      
!     A trapezoid shape opening calculation
!     XX1: Start timing of suction
!     XX2: Till this timing the opening is fully opened
!     XX3: Fully opening is till this timing
!     XX4: The end of suction
!     XX5: Start timing of discharge
!     XX6: Same with suction ...
      if ((xt.ge.0.0d0) .and. (xt.le.xx(4))) then
         if (xt.le.xx(2)) then
            XK1=(XT-XX(1))/(XX(2)-XX(1))
            B=pi*xk1
         else if (xt.le.xx(3)) then
            B=pi
         else if (xt.le.xx(4)) then
            XK1=(XT-XX(3))/(XX(4)-XX(3))
            B=pi-pi*xk1
         end if
      end if
      if (xt.gt.xx(4)) then
         b = 0.0d0
      end if
      cv(1) = cv10*B

      if ((xt.ge.xx(5)).and.(xt.le.xx(8))) then
         if (xt.le.xx(6)) then
            XK2=(XT-XX(5))/(XX(6)-XX(5))
            B=pi*xk2
         else if (xt.le.xx(7)) then
            B=pi
         else if (xt.le.xx(8)) then
            XK2=(XT-XX(7))/(XX(8)-XX(7))
            B=pi-pi*xk2
         end if
      end if
      if (xt.lt.xx(5)) then
         b = 0.0d0
      end if
      cv(2) = cv20*B

      end subroutine
