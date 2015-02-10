c-----
c SBML Model : potassium                                              
c              ~~~~~~~~~
c       Date : Wed Jan 21 19:15:19 2015
c              
c     Author : automated transcription by 'sbml2fortran'
c              
c Copyright (C) Zuse Institute Berlin, CSB Group
c-----
c
      subroutine ydot_LIMEX ( n, nz, t, y_state, dy, b, ir, ic, info )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer             j, n, nz, info
c
      integer             ir(*), ic(*)
c
      double precision    t, y_state(n), dy(n), b(*), one, zero
c
      parameter           ( one = 1.0D0, zero = 0.0D0 )
c
c-----------------------------------------------------------------------
c
      nz = n
c
      do j = 1, nz
c
         ir(j) = j
         ic(j) = j
         b(j) = one
c
      end do
c
c 
c
c-----------------------------------------------------------------------
c
      call ydot ( n, t, y_state, dy, info )
c
c-----------------------------------------------------------------------
c
      return
      end subroutine ydot_LIMEX
c
c=======================================================================
c
      subroutine ydot ( n, t, y_state, dy, info )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer             j, n, info
c
      double precision    t, y_state(n), dy(n), zero
c
c     ---
c     com(  1)  :         default                                   
c     com(  2)  :         c1                                        
c     com(  3)  :         c2                                        
c     ---
c
c     ---
c     spe(  1)  :         K_ECF                                     
c     spe(  2)  :         K_ICF                                     
c     spe(  3)  :         K_urin                                    
c     spe(  4)  :         K_git                                     
c     spe(  5)  :         K_tiss                                    
c     spe(  6)  :         K_sal                                     
c     spe(  7)  :         K_feed                                    
c     spe( 11)  :         Insulin                                   
c     spe( 12)  :         s20                                       
c     spe( 13)  :         s23                                       
c     spe( 14)  :         s27                                       
c     spe( 16)  :         s29                                       
c     spe( 22)  :         s35                                       
c     ---
c
c     ---
c     par(  1)  :         global_p1                                 
c     par(  2)  :         global_p2                                 
c     par(  3)  :         global_p3                                 
c     par(  4)  :         global_p4                                 
c     par(  5)  :         global_p5                                 
c     par(  6)  :         global_p6                                 
c     par(  7)  :         global_p7                                 
c     par(  8)  :         global_p8                                 
c     par(  9)  :         global_p9                                 
c     par( 10)  :         global_p10                                
c     par( 11)  :         global_p11                                
c     par( 12)  :         global_p12                                
c     par( 13)  :         global_p13                                
c     par( 14)  :         global_p14                                
c     par( 15)  :         global_p15                                
c     par( 16)  :         global_p16                                
c     par( 17)  :         global_p17                                
c     par( 18)  :         global_p18                                
c     par( 19)  :         global_p19                                
c     par( 20)  :         global_p20                                
c     par( 21)  :         global_p21                                
c     par( 22)  :         global_p22                                
c     par( 23)  :         global_p23                                
c     par( 24)  :         global_p24                                
c     par( 25)  :         global_p25                                
c     par( 26)  :         global_p26                                
c     par( 27)  :         global_p27                                
c     par( 28)  :         global_p28                                
c     par( 29)  :         global_p29                                
c     par( 30)  :         global_p30                                
c     par( 31)  :         global_p31                                
c     par( 32)  :         global_p32                                
c     par( 33)  :         global_p33                                
c     par( 34)  :         global_p34                                
c     par( 35)  :         global_p35                                
c     par( 36)  :         global_pH                                 
c     par( 37)  :         global_p39                                
c     par( 38)  :         global_p40                                
c     par( 39)  :         global_p41                                
c     par( 40)  :         global_p42                                
c     par( 41)  :         global_p43                                
c     par( 42)  :         global_p44                                
c     par( 43)  :         global_p45                                
c     par( 44)  :         global_p46                                
c     par( 45)  :         global_p47                                
c     par( 46)  :         global_p48                                
c     par( 47)  :         global_p49                                
c     par( 48)  :         global_p50                                
c     ---
c
      double precision    com(3)
      double precision    spe(22)
      double precision    par(48)
      double precision    rul(1)
      double precision    rea(20)
      logical             eve(0)
c
      common              /SBMLVARIABLES/ com, spe, par, rul, rea, eve
c
c
c 
c 
      double precision    abs
      double precision    acos, arccos, acosh, arccosh
      double precision    acot, arccot, acoth, arccoth
      double precision    acsc, arccsc, acsch, arccsch
      double precision    asec, arcsec, asech, arcsech
      double precision    asin, arcsin, asinh, arcsinh
      double precision    atan, arctan, atanh, arctanh
      double precision    ceil, ceiling, floor, factorial
      double precision    cos, cosh, cot, coth, csc, csch
      double precision    exp, ln, log, log10
      double precision    pow, power, root, sqr, sqrt
      double precision    sec, sech, sin, sinh, tan, tanh
      logical             and, or, not
      logical             eq, neq, gt, geq, lt, leq
c
      parameter           ( zero = 0.0D0 )
c
c-----------------------------------------------------------------------
c
      info = 0
c
c-----------------------------------------------------------------------
c
      do j = 1, n
c
         spe(j) = y_state(j)
c
         dy(j) = zero
c
      end do
c
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      rul(1) = 7.500000D+00 - spe(13) / 4.000000D+01
c 
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      rea(1) = (1.000000D+00 + 1.000000D+00 / (1.000000D+00 +
     &          pow(spe(1) / par(29), 1.000000D+01))) * par(5) *
     &          (1.000000D+00 + 1.000000D+00 / (1.000000D+00 +
     &          pow(rul(1) / par(20), 1.000000D+01))) * pow(spe(2),
     &          2.000000D+00) / (pow(spe(2), 2.000000D+00) +
     &          pow(par(23), 2.000000D+00)) * (1.000000D+00 + par(19)
     &          * pow(spe(2), 2.000000D+00) / (pow(spe(2),
     &          2.000000D+00) + pow(par(14), 2.000000D+00))) /
     &          par(34)
      rea(2) = (par(8) + par(9) * pow(spe(11), 1.000000D+01) /
     &          (pow(spe(11), 1.000000D+01) + pow(par(3),
     &          1.000000D+01))) * pow(spe(1), 2.000000D+00) /
     &          (pow(spe(1), 2.000000D+00) + pow(par(33),
     &          2.000000D+00)) * (1.000000D+00 + par(21) *
     &          pow(rul(1), 1.000000D+01) / (pow(rul(1),
     &          1.000000D+01) + pow(par(20), 1.000000D+01))) /
     &          par(34)
      rea(3) = par(18) * spe(1) / par(34)
      rea(4) = par(4) * spe(1) / par(34)
      rea(5) = pow(spe(4), 1.000000D+01) / (pow(spe(4),
     &          1.000000D+01) + pow(par(15), 1.000000D+01)) * par(25)
     &          * spe(1) / par(34)
      rea(6) = 1.000000D+00 / (1.000000D+00 + pow(spe(4) / par(15),
     &          1.000000D+01)) * par(26) * (par(27) - spe(1)) *
     &          pow(spe(5), 1.000000D+01) / (pow(spe(5),
     &          1.000000D+01) + pow(1.468375D+03, 1.000000D+01)) /
     &          par(34)
      rea(7) = (1.000000D+00 + par(13) * pow(spe(1), 5.000000D+00) /
     &          (pow(spe(1), 5.000000D+00) + pow(par(24),
     &          5.000000D+00))) * par(6) * spe(4) * (1.000000D+00 +
     &          par(16) * spe(1) * pow(spe(1), 1.000000D+01) /
     &          (pow(spe(1), 1.000000D+01) + pow(par(22),
     &          1.000000D+01)))
      rea(8) = par(31) * spe(4)
      rea(9) = par(32) * spe(6)
      rea(10) = par(48) * spe(11)
      rea(11) = par(3) + par(2) * spe(16)
      rea(12) = par(30) * spe(7)
      rea(13) = par(43) + par(40) * pow(spe(12), 2.000000D+00) /
     &          (pow(spe(12), 2.000000D+00) + pow(par(41),
     &          2.000000D+00))
      rea(14) = par(42) * spe(13)
      rea(15) = par(37) * spe(12)
      rea(16) = par(44) * spe(12)
      rea(17) = 3.000000D+00 + par(46) * spe(14)
      rea(18) = par(45) * spe(16)
      rea(19) = par(47) * spe(14)
      rea(20) = 4.000000D+01 * 3.141600D+00 * cos((t - 9.000000D+00)
     &          * 3.141600D+00 / 1.200000D+01)
c 
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      dy(1) = ( + rea(1) - rea(2) - rea(3) - rea(4) - rea(5) +
     &          rea(6) - rea(7) + rea(8) ) / com(3)
      dy(2) = ( - rea(1) + rea(2) ) / com(2)
      dy(3) = ( + rea(7) ) / com(1)
      dy(4) = ( - rea(8) + rea(9) + rea(12) ) / com(1)
      dy(5) = ( + rea(5) - rea(6) ) / com(1)
      dy(6) = ( + rea(3) - rea(9) ) / com(1)
      dy(7) = ( - rea(12) + rea(15) ) / com(1)
      dy(11) = ( - rea(10) + rea(11) ) / com(3)
      dy(12) = ( + rea(20) ) / com(1)
      dy(13) = ( + rea(13) - rea(14) ) / com(1)
      dy(14) = ( + rea(16) - rea(17) - rea(19) ) / com(1)
      dy(16) = ( + rea(17) - rea(18) ) / com(3)
      dy(22) = ( - rea(20) ) / com(1)
c 
c-----------------------------------------------------------------------
c
      return
      end subroutine ydot
c
c=======================================================================
c
      subroutine check_events ( t, n, ev )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    t
c
      integer             n
c
      logical             ev(n)
c
      double precision    com(3)
      double precision    spe(22)
      double precision    par(48)
      double precision    rul(1)
      double precision    rea(20)
      logical             eve(0)
c
      common              /SBMLVARIABLES/ com, spe, par, rul, rea, eve
c
c
c 
c 
      double precision    abs
      double precision    acos, arccos, acosh, arccosh
      double precision    acot, arccot, acoth, arccoth
      double precision    acsc, arccsc, acsch, arccsch
      double precision    asec, arcsec, asech, arcsech
      double precision    asin, arcsin, asinh, arcsinh
      double precision    atan, arctan, atanh, arctanh
      double precision    ceil, ceiling, floor, factorial
      double precision    cos, cosh, cot, coth, csc, csch
      double precision    exp, ln, log, log10
      double precision    pow, power, root, sqr, sqrt
      double precision    sec, sech, sin, sinh, tan, tanh
      logical             and, or, not
      logical             eq, neq, gt, geq, lt, leq
c
c-----------------------------------------------------------------------
c
c 
c-----------------------------------------------------------------------
c
      return
      end subroutine check_events
c
c=======================================================================
c
      subroutine init_ode ( c, ic, lc ,  s, is, ls ,  p, ip, lp )
c
      implicit none
c
c-----------------------------------------------------------------------
c
c     ---
c     com(  1)  :         default                                   
c     com(  2)  :         c1                                        
c     com(  3)  :         c2                                        
c     ---
c
c     ---
c     spe(  1)  :         K_ECF                                     
c     spe(  2)  :         K_ICF                                     
c     spe(  3)  :         K_urin                                    
c     spe(  4)  :         K_git                                     
c     spe(  5)  :         K_tiss                                    
c     spe(  6)  :         K_sal                                     
c     spe(  7)  :         K_feed                                    
c     spe( 11)  :         Insulin                                   
c     spe( 12)  :         s20                                       
c     spe( 13)  :         s23                                       
c     spe( 14)  :         s27                                       
c     spe( 16)  :         s29                                       
c     spe( 22)  :         s35                                       
c     ---
c
c     ---
c     par(  1)  :         global_p1                                 
c     par(  2)  :         global_p2                                 
c     par(  3)  :         global_p3                                 
c     par(  4)  :         global_p4                                 
c     par(  5)  :         global_p5                                 
c     par(  6)  :         global_p6                                 
c     par(  7)  :         global_p7                                 
c     par(  8)  :         global_p8                                 
c     par(  9)  :         global_p9                                 
c     par( 10)  :         global_p10                                
c     par( 11)  :         global_p11                                
c     par( 12)  :         global_p12                                
c     par( 13)  :         global_p13                                
c     par( 14)  :         global_p14                                
c     par( 15)  :         global_p15                                
c     par( 16)  :         global_p16                                
c     par( 17)  :         global_p17                                
c     par( 18)  :         global_p18                                
c     par( 19)  :         global_p19                                
c     par( 20)  :         global_p20                                
c     par( 21)  :         global_p21                                
c     par( 22)  :         global_p22                                
c     par( 23)  :         global_p23                                
c     par( 24)  :         global_p24                                
c     par( 25)  :         global_p25                                
c     par( 26)  :         global_p26                                
c     par( 27)  :         global_p27                                
c     par( 28)  :         global_p28                                
c     par( 29)  :         global_p29                                
c     par( 30)  :         global_p30                                
c     par( 31)  :         global_p31                                
c     par( 32)  :         global_p32                                
c     par( 33)  :         global_p33                                
c     par( 34)  :         global_p34                                
c     par( 35)  :         global_p35                                
c     par( 36)  :         global_pH                                 
c     par( 37)  :         global_p39                                
c     par( 38)  :         global_p40                                
c     par( 39)  :         global_p41                                
c     par( 40)  :         global_p42                                
c     par( 41)  :         global_p43                                
c     par( 42)  :         global_p44                                
c     par( 43)  :         global_p45                                
c     par( 44)  :         global_p46                                
c     par( 45)  :         global_p47                                
c     par( 46)  :         global_p48                                
c     par( 47)  :         global_p49                                
c     par( 48)  :         global_p50                                
c     ---
c
c-----------------------------------------------------------------------
c
      double precision    com(3)
      double precision    spe(22)
      double precision    par(48)
      double precision    rul(1)
      double precision    rea(20)
      logical             eve(0)
c
      common              /SBMLVARIABLES/ com, spe, par, rul, rea, eve
c
      double precision    c(*), s(*), p(*)
c
      integer             ic(*), is(*), ip(*), lc, ls, lp
c
      integer             mc, ms, mp, j
c
c-----------------------------------------------------------------------
c
      com(  1) = 1.000000D+00
      com(  2) = 1.993800D+00
      com(  3) = 8.916000D-01
c 
c
      spe(  1) = 3.430000D+00
      spe(  2) = 2.542000D+01
      spe(  3) = 0.000000D+00
      spe(  4) = 3.073000D+00
      spe(  5) = 1.509600D+03
      spe(  6) = 4.340000D-01
      spe(  7) = 0.000000D+00
      spe( 11) = 2.200000D+01
      spe( 12) = 1.405887D+02
      spe( 13) = 3.000000D+00
      spe( 14) = 6.250000D+01
      spe( 16) = 3.730000D+00
      spe( 22) = 0.000000D+00
c 
c
      par(  1) = 3.730000D+00
      par(  2) = 3.400000D+01
      par(  3) = 2.450000D+01
      par(  4) = 2.000000D-04
      par(  5) = 1.404000D+01
      par(  6) = 5.000000D-02
      par(  7) = 7.500000D+01
      par(  8) = 1.491000D+01
      par(  9) = 9.990000D+00
      par( 10) = 6.000000D+02
      par( 11) = 3.430000D+00
      par( 12) = 2.542000D+01
      par( 13) = 6.010000D+00
      par( 14) = 4.000000D+00
      par( 15) = 2.200000D+01
      par( 16) = 3.000000D+00
      par( 17) = 3.090000D-01
      par( 18) = 8.000000D-02
      par( 19) = 1.373000D-01
      par( 20) = 7.380000D+00
      par( 21) = 1.085000D-01
      par( 22) = 8.000000D+00
      par( 23) = 1.500000D+01
      par( 24) = 5.000000D+00
      par( 25) = 2.898000D-01
      par( 26) = 2.150000D+00
      par( 27) = 5.950000D+00
      par( 28) = 1.100000D-03
      par( 29) = 3.000000D+00
      par( 30) = 9.500000D-01
      par( 31) = 3.500000D-01
      par( 32) = 8.000000D-01
      par( 33) = 2.000000D+00
      par( 34) = 1.121600D+00
      par( 35) = 2.236240D+00
      par( 36) = 0.000000D+00
      par( 37) = 1.135000D-02
      par( 38) = 0.000000D+00
      par( 39) = 7.312500D+02
      par( 40) = 4.500000D+00
      par( 41) = 1.000000D+01
      par( 42) = 1.200000D+00
      par( 43) = 3.000000D+00
      par( 44) = 8.540000D-02
      par( 45) = 1.000000D+00
      par( 46) = 1.500000D-02
      par( 47) = 9.850000D-01
      par( 48) = 7.000000D+00
c 
c
c 
c
c-----------------------------------------------------------------------
c
      mc = min0(3,lc)
c
      if ( ic(1).eq.0 ) then
c
         do j = 1, mc
c
            com( j ) = c(j)
c
         end do
c
      else
c
         do j = 1, mc
c
            com( ic(j) ) = c(j)
c
         end do
c
      end if
c
c
      if ( lc.eq.0 ) then
c
         lc = 3
c
      end if
c
c-----------------------------------------------------------------------
c
      ms = min0(22,ls)
c
      if ( is(1).eq.0 ) then
c
         do j = 1, ms
c
            spe( j ) = s(j)
c
         end do
c
      else
c
         do j = 1, ms
c
            spe( is(j) ) = s(j)
c
         end do
c
      end if
c
c
      if ( ls.eq.0 ) then
c
         ls = 22
c
      end if
c
c-----------------------------------------------------------------------
c
      mp = min0(48,lp)
c
      if ( ip(1).eq.0 ) then
c
         do j = 1, mp
c
            par( j ) = p(j)
c
         end do
c
      else
c
         do j = 1, mp
c
            par( ip(j) ) = p(j)
c
         end do
c
      end if
c
c
      if ( lp.eq.0 ) then
c
         lp = 48
c
      end if
c
c-----------------------------------------------------------------------
c
      return
      end subroutine init_ode
c
c=======================================================================
c
      subroutine get_compartment_ids ( idc, ncom )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character           idc(3)*(*)
c
      integer             ncom
c
c-----------------------------------------------------------------------
c
      ncom = 3
c     ---
      idc(  1) = 'default'//char(0)
      idc(  2) = 'c1'//char(0)
      idc(  3) = 'c2'//char(0)
c     ---
c
c-----------------------------------------------------------------------
c
      return
      end subroutine get_compartment_ids
c
c=======================================================================
c
      subroutine get_species_ids ( ids, nspe )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character           ids(22)*(*)
c
      integer             j, nspe
c
c-----------------------------------------------------------------------
c
      nspe = 22
c     ---
      do j = 1, nspe
c
        ids(j) = char(0)
c
      end do
c     ---
      ids(  1) = 'K_ECF'//char(0)
      ids(  2) = 'K_ICF'//char(0)
      ids(  3) = 'K_urin'//char(0)
      ids(  4) = 'K_git'//char(0)
      ids(  5) = 'K_tiss'//char(0)
      ids(  6) = 'K_sal'//char(0)
      ids(  7) = 'K_feed'//char(0)
      ids( 11) = 'Insulin'//char(0)
      ids( 12) = 's20'//char(0)
      ids( 13) = 's23'//char(0)
      ids( 14) = 's27'//char(0)
      ids( 16) = 's29'//char(0)
      ids( 22) = 's35'//char(0)
c     ---
c
c-----------------------------------------------------------------------
c
      return
      end subroutine get_species_ids
c
c=======================================================================
c
      subroutine get_parameter_ids ( idp, npar )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character           idp(48)*(*)
c
      integer             j, npar
c
c-----------------------------------------------------------------------
c
      npar = 48
c     ---
      do j = 1, npar
c
        idp(j) = char(0)
c
      end do
c     ---
      idp(  1) = 'global_p1'//char(0)
      idp(  2) = 'global_p2'//char(0)
      idp(  3) = 'global_p3'//char(0)
      idp(  4) = 'global_p4'//char(0)
      idp(  5) = 'global_p5'//char(0)
      idp(  6) = 'global_p6'//char(0)
      idp(  7) = 'global_p7'//char(0)
      idp(  8) = 'global_p8'//char(0)
      idp(  9) = 'global_p9'//char(0)
      idp( 10) = 'global_p10'//char(0)
      idp( 11) = 'global_p11'//char(0)
      idp( 12) = 'global_p12'//char(0)
      idp( 13) = 'global_p13'//char(0)
      idp( 14) = 'global_p14'//char(0)
      idp( 15) = 'global_p15'//char(0)
      idp( 16) = 'global_p16'//char(0)
      idp( 17) = 'global_p17'//char(0)
      idp( 18) = 'global_p18'//char(0)
      idp( 19) = 'global_p19'//char(0)
      idp( 20) = 'global_p20'//char(0)
      idp( 21) = 'global_p21'//char(0)
      idp( 22) = 'global_p22'//char(0)
      idp( 23) = 'global_p23'//char(0)
      idp( 24) = 'global_p24'//char(0)
      idp( 25) = 'global_p25'//char(0)
      idp( 26) = 'global_p26'//char(0)
      idp( 27) = 'global_p27'//char(0)
      idp( 28) = 'global_p28'//char(0)
      idp( 29) = 'global_p29'//char(0)
      idp( 30) = 'global_p30'//char(0)
      idp( 31) = 'global_p31'//char(0)
      idp( 32) = 'global_p32'//char(0)
      idp( 33) = 'global_p33'//char(0)
      idp( 34) = 'global_p34'//char(0)
      idp( 35) = 'global_p35'//char(0)
      idp( 36) = 'global_pH'//char(0)
      idp( 37) = 'global_p39'//char(0)
      idp( 38) = 'global_p40'//char(0)
      idp( 39) = 'global_p41'//char(0)
      idp( 40) = 'global_p42'//char(0)
      idp( 41) = 'global_p43'//char(0)
      idp( 42) = 'global_p44'//char(0)
      idp( 43) = 'global_p45'//char(0)
      idp( 44) = 'global_p46'//char(0)
      idp( 45) = 'global_p47'//char(0)
      idp( 46) = 'global_p48'//char(0)
      idp( 47) = 'global_p49'//char(0)
      idp( 48) = 'global_p50'//char(0)
c     ---
c
c-----------------------------------------------------------------------
c
      return
      end subroutine get_parameter_ids
c
c=======================================================================
c
      subroutine get_model_ids ( idm, nid )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      character           idm(4)*(*)
c
      integer             j, nid
c
c-----------------------------------------------------------------------
c
      nid = 4
c     ---
      do j = 1, nid
c
        idm(j) = char(0)
c
      end do
c     ---
      idm(  1) = 'potassium'//char(0)
      idm(  2) = 'Wed Jan 21 19:15:19 2015'//char(0)
      idm(  3) = '001421864119'//char(0)
      idm(  4) = 'potmod_2_sin_dmi.xml'//char(0)
c     ---
c
c-----------------------------------------------------------------------
c
      return
      end subroutine get_model_ids
c
c=======================================================================
c
c 
c
c 
c

