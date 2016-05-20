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

!     Valve timing information is needed here because otherwise
!     the opening data being written below will be strange.
      call valve(x, prmt)

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
