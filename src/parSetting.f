      subroutine parSetting(y, dery, prmt, ndim)
!     Global variables used by other subroutines are calculated here.
      use data
      implicit none

!     These variables are used by the Runge-Kutta subroutine.
      double precision y, dery, prmt
      DIMENSION Y(3),DERY(3),PRMT(8)
      integer ndim

      CV(1)=B1*21.2
      CV(2)=B2*21.2
      CV(3)=B3
      PRMT(1)=0.0
      PRMT(4)=0.0001
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
!     R=2.077,(for He [MPa,cm3/g/K])
      Y(1)=0.0
      Y(2)=0.0
      Y(3)=0.0
      DERY(1)=0.3
      DERY(2)=0.3
      DERY(3)=0.4
      NDIM=3

      end subroutine
