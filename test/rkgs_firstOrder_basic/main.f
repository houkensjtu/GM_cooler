      program main
!     A program used to test the ability of the fortran Runge-Kutta subroutine.
!     Only one differential equation is solved here, y'(x) = x
!     Thus the aux, y, dery has only dimension 1.
!     Find http://media.ibm1130.org/1130-006-ocr.pdf for further explaination.
      implicit none
      double precision prmt, y, dery, aux, x
      dimension prmt(8), aux(8,1), y(1), dery(1)
      
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
      prmt(2) = 10.0d0
      prmt(3) = 0.01d0
      prmt(4) = 0.01d0
      
!     Should not forget to set the initial value of y.
      y(1) = 0.0d0

!     Ony one equation here, so ndim=1.
      ndim = 1
      ihlf = 3

!     Initial value of dery(i) is used as a weight factor for each equation.
!     While this value will be replaced by the value of y'(x) during rkgs.
      dery(1) = 1.0d0
      
      call rkgs(prmt, y, dery, ndim, ihlf, aux)
      end program

      
!     User should provide the fct and outp subroutine by themselves.
      subroutine fct(x, y, dery)
      implicit none

!     The preferred way to declare an array is to first declare its type,
!     then following by its dimension.
      double precision x, y ,dery
      dimension dery(1),y(1)
      
!     Calculate the y'(x) as you want.
      dery(1) =  x
      end subroutine
      
      
      subroutine outp(x, y, dery, ihlf, ndim, prmt)
      implicit none
      double precision prmt, y, dery, aux, x
      dimension prmt(8), aux(8,1), dery(1), y(1)
      integer ihlf, ndim

!     Print the result as you want.
      print *,x,y(1),dery(1)

      end subroutine

