      function pressureAssist(x,f)
      use data, only: PH, PL
      implicit none
      
      double precision pressureAssist
      double precision x, f

      if ( (mod(x, 1.0D0/f).gt.(0.2D0/f)) .and. 
     + (mod(x,1.0d0/f).lt.(0.85d0/f))) then
          pressureAssist = PL
      else 
          pressureAssist = PH
      end if

      end function
