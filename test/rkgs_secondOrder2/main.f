      program main
!     A program used to test the ability of the fortran Runge-Kutta subroutine.
!     A second-order differential equation is solved here, y''(x) - 2y'(x) + 10y(x) = 0.
!     Thus the aux, y, dery has the dimension of 2.

!     The second order equation is solved by splitting it into 2 first-order:
!     y1'(t) = y2
!     y2'(t) = 2*y2 - 10*y1
!     So y1''(t) - 2y1'(t) + 10y1(t) = 0
!     The exact solution of y1(t) should be y1(t) = e^x*(4cos3x - sin3x)

!     Find http://media.ibm1130.org/1130-006-ocr.pdf for further explaination.

      implicit none
      double precision prmt, y, dery, aux, x
      dimension prmt(8), aux(8,2), y(2), dery(2)
      
!     ihlf is a auxilary value, canbe omitted.
!     ndim is the number of equations.
      integer ihlf, ndim

!     prmt(1) is the begining value of x.
!     prmt(2) is the end value of x.
!     prmt(3) is the initial increment.
!     prmt(4) is the upper error bound.
!     The rest of prmt could be ommited.
!     Add d0 following the number to avoid double precision error.
      prmt(1) = 0.0d0
      prmt(2) = 1.0d0
      prmt(3) = 0.1d0
      prmt(4) = 0.0001d0
      
!     Should not forget to set the initial value of y.
      y(1) = 4.0d0
      y(2) = 1.0d0

!     Two equations here, so ndim=2.
      ndim = 2
      ihlf = 3

!     Initial value of dery(i) is used as a weight factor for each equation.
!     While this value will be replaced by the value of y'(x) during rkgs.
      dery(1) = 0.5d0
      dery(2) = 0.5d0
      
      call rkgs(prmt, y, dery, ndim, ihlf, aux)
      end program

      
!     User should provide the fct and outp subroutine by themselves.
      subroutine fct(x, y, dery)
      implicit none

!     The preferred way to declare an array is to first declare its type,
!     then following by its dimension.
      double precision x, y ,dery
      dimension dery(2),y(2)
      
!     Calculate the y'(x) as you want.
      dery(1) = y(2)
      dery(2) = 2.0D0*y(2) - 10.0d0*y(1)
      end subroutine
      
      
      subroutine outp(x, y, dery, ihlf, ndim, prmt)
      implicit none
      double precision prmt, y, dery, aux, x
      dimension prmt(8), aux(8,2), dery(2), y(2)
      integer ihlf, ndim
      double precision y_sol

!     Print the result as you want.
      if (x.eq.0.0D0) then
         write(*,'(6A10)') 'x','y1','y2','dery1','dery2','y_sol'
      end if

!     Calculate the exact solution.
      y_sol = 2.718**x * ( 4.0d0*cos(3.0d0*x) - sin(3.0d0*x) )

!     Output
      write(*,100) x,y(1),y(2),dery(1),dery(2),y_sol
 100  format(6f10.3)
      end subroutine
