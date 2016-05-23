      function pressureAssist(x,f)
      use data, only: PH, PL
      implicit none
      
      double precision pressureAssist
      double precision x, f

      if (mod(x, 1.0D0/f).gt.(0.5D0/f)) then
          pressureAssist = PH
      else 
          pressureAssist = PL
      end if

      end function
