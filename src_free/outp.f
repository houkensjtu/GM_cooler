      SUBROUTINE OUTP(X,Y,DERY,IHLF,NDIM,PRMT)
!     This subroutine is a "MUST" for Runge-Kutta subroutine.
!     Look http://media.ibm1130.org/1130-006-ocr.pdf for detailed description.
!     However, instead of writing out data, this subroutine is mainly
!     aimed at updating the value of Valve opening - CV(3).

      use data
      implicit none
      double precision y, dery, prmt, aux, x
      DIMENSION Y(5),DERY(5),AUX(8,5),PRMT(8)
      integer ndim, ihlf

!     Valve timing information is needed here because otherwise
!     the opening data being written below will be strange.

!     valve() is the original version (Triangle opening)
!     valveTrapezoid is the new trapezoid shape opening
!     call valve(x, prmt)
      call valveTrapezoid(x, prmt)

      if (output_mode == '1') then
         output_timing = prmt(7)
      else if (output_mode == '2') then
         output_timing = 0.0d0
      end if

      IF(x.ge.output_timing) then

         log_step = log_step + 1

!     Expansion room volume
         D(1,log_step)=V2

!     Compression room pressure.
         D(2,log_step)=P1

!     Expansion room pressure.
         D(3,log_step)=P2

!     mass flow rate at regenerator
         D(4,log_step)=DERY(3)

!     mass flow rate at Valve1
         D(5,log_step)=DERY(1)

!     mass flow rate at Valve2
         D(6,log_step)=DERY(2)

!     Valve 1 opening
         D(7,log_step)=CV(1)/21.2

!     Valve 2 opening
         D(8,log_step)=CV(2)/21.2

!     Place holder for displacer calculation.
!     pa is the pressure of assist room.
         D(9,log_step)=pa

!     Position of displacer.
         D(10,log_step)=y(4)

!     Velocity of displacer.
         D(11,log_step)=dery(4)

      ENDIF

      RETURN
      END SUBROUTINE
