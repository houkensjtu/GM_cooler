      program main

      double precision prmt, y, dery, aux, x
      dimension prmt(8), aux(8)
      integer ihlf, ndim
          print *,ihlf

      prmt(1) = 0.0
      prmt(2) = 2.0
      prmt(3) = 0.01
      prmt(4) = 0.001
      
      ndim = 1
      ihlf = 3

      call rkgs(prmt, y, dery, ndim, ihlf, aux)
      end program

      subroutine fct(x, y, dery)
      
      double precision x, y ,dery
      dery = x
      end subroutine
      
      
      subroutine outp(x, y, dery, ihlf, ndim, prmt)

      double precision prmt, y, dery, aux, x
      dimension prmt(8), aux(8)
      integer ihlf, ndim
      
      print *,aux
      end subroutine

