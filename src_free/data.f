      module data
!     Variables declaration seperated from source code.
!     Instead of old style "IMPLICIT" definition, all variables
!     are forced to be defined before using.
      implicit none
      double precision TH,TRH,TRE,TE,R,P1,P2,PH,PL,XX,CV,D
      double precision VT,VRH,VRE,C1,C2,CV10,CV20,PHAI,V1,V2,OMEGA
      double precision XT
      
      double precision vrd, f, fie, b1, b2, regFlowCoef

!     8 timing parameter used to form a Trapezoid shape opening
      double precision x1, x2, x3, x4, x5, x6, x7, x8
      double precision phai0, pi, cp, tr, pm
      double precision v20, v10
      
      double precision amr, ame, fm1, fm2, wf0, wf1, wu1, wf2, wf10
      double precision wf20, we1, we2, q, hr, wact, wee, wu2
      double precision caract, carnot, cmax, cmax2, dp

      double precision pp1, pp2, ppx1, ppx2, p0
      double precision ff, xk1, xk2, xk3, xk4, xk5, xk6, xk7, xk8
      double precision c, s, b

      integer j, mode, jd, jj, k, i
      integer nx, nf, n
      integer log_step
      double precision k_float
      dimension d(12,100000), xx(8), cv(3)
      double precision am0, amemax

!     Assist pressure.
      double precision pa
      
      character(len=3) :: drive_mode, output_mode
      double precision output_timing

      double precision displacerMass, dampCoef, springCoef
      end module data
