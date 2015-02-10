c-----
c SBML Model : PAEON_V2                                               
c              ~~~~~~~~
c       Date : Mon Feb  2 18:21:20 2015
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
c     ---
c
c     ---
c     spe(  1)  :         LH_pit                                    
c     spe(  2)  :         LH_blood                                  
c     spe(  3)  :         R_LH                                      
c     spe(  4)  :         LH_R                                      
c     spe(  5)  :         R_LH_des                                  
c     spe(  6)  :         FSH_pit                                   
c     spe(  7)  :         FSH_blood                                 
c     spe(  8)  :         R_FSH                                     
c     spe(  9)  :         FSH_R                                     
c     spe( 10)  :         R_FSH_des                                 
c     spe( 11)  :         S_foll                                    
c     spe( 12)  :         AF1                                       
c     spe( 13)  :         AF2                                       
c     spe( 14)  :         AF3                                       
c     spe( 15)  :         AF4                                       
c     spe( 16)  :         PrF                                       
c     spe( 17)  :         OvF                                       
c     spe( 18)  :         Sc1                                       
c     spe( 19)  :         Sc2                                       
c     spe( 20)  :         Lut1                                      
c     spe( 21)  :         Lut2                                      
c     spe( 22)  :         Lut3                                      
c     spe( 23)  :         Lut4                                      
c     spe( 24)  :         E2                                        
c     spe( 25)  :         P4                                        
c     spe( 26)  :         InhA                                      
c     spe( 27)  :         InhB                                      
c     spe( 28)  :         InhA_tau                                  
c     spe( 29)  :         GnRH_G                                    
c     spe( 30)  :         R_GnRH_a                                  
c     spe( 31)  :         R_GnRH_i                                  
c     spe( 32)  :         GnRH_R_a                                  
c     spe( 33)  :         GnRH_R_i                                  
c     ---
c
c     ---
c     par(  1)  :         global_p_001_001                          
c     par(  2)  :         global_p_001_002                          
c     par(  3)  :         global_p_001_003                          
c     par(  4)  :         global_p_001_004                          
c     par(  5)  :         global_p_001_005                          
c     par(  6)  :         global_p_001_006                          
c     par(  7)  :         global_p_001_007                          
c     par(  8)  :         global_p_001_008                          
c     par(  9)  :         global_p_001_009                          
c     par( 10)  :         global_p_001_010                          
c     par( 11)  :         global_p_002_001                          
c     par( 12)  :         global_p_002_002                          
c     par( 13)  :         global_p_002_003                          
c     par( 14)  :         global_p_003_001                          
c     par( 15)  :         global_p_004_001                          
c     par( 16)  :         global_p_006_001                          
c     par( 17)  :         global_p_006_002                          
c     par( 18)  :         global_p_006_003                          
c     par( 19)  :         global_p_006_004                          
c     par( 20)  :         global_p_006_005                          
c     par( 21)  :         global_p_006_006                          
c     par( 22)  :         global_p_006_007                          
c     par( 23)  :         global_p_006_008                          
c     par( 24)  :         global_p_006_009                          
c     par( 25)  :         global_p_006_010                          
c     par( 26)  :         global_p_006_011                          
c     par( 27)  :         global_p_007_001                          
c     par( 28)  :         global_p_007_002                          
c     par( 29)  :         global_p_008_001                          
c     par( 30)  :         global_p_009_001                          
c     par( 31)  :         global_p_011_001                          
c     par( 32)  :         global_p_011_002                          
c     par( 33)  :         global_p_011_003                          
c     par( 34)  :         global_p_011_004                          
c     par( 35)  :         global_p_011_005                          
c     par( 36)  :         global_p_011_006                          
c     par( 37)  :         global_p_011_007                          
c     par( 38)  :         global_p_012_001                          
c     par( 39)  :         global_p_012_002                          
c     par( 40)  :         global_p_012_003                          
c     par( 41)  :         global_p_012_004                          
c     par( 42)  :         global_p_013_001                          
c     par( 43)  :         global_p_013_002                          
c     par( 44)  :         global_p_013_003                          
c     par( 45)  :         global_p_014_001                          
c     par( 46)  :         global_p_014_002                          
c     par( 47)  :         global_p_014_003                          
c     par( 48)  :         global_p_014_004                          
c     par( 49)  :         global_p_015_001                          
c     par( 50)  :         global_p_015_002                          
c     par( 51)  :         global_p_015_003                          
c     par( 52)  :         global_p_016_001                          
c     par( 53)  :         global_p_016_002                          
c     par( 54)  :         global_p_017_001                          
c     par( 55)  :         global_p_017_002                          
c     par( 56)  :         global_p_017_003                          
c     par( 57)  :         global_p_017_004                          
c     par( 58)  :         global_p_018_001                          
c     par( 59)  :         global_p_018_002                          
c     par( 60)  :         global_p_018_003                          
c     par( 61)  :         global_p_018_004                          
c     par( 62)  :         global_p_019_001                          
c     par( 63)  :         global_p_020_001                          
c     par( 64)  :         global_p_020_002                          
c     par( 65)  :         global_p_020_003                          
c     par( 66)  :         global_p_020_004                          
c     par( 67)  :         global_p_021_001                          
c     par( 68)  :         global_p_022_001                          
c     par( 69)  :         global_p_023_001                          
c     par( 70)  :         global_p_024_001                          
c     par( 71)  :         global_p_024_002                          
c     par( 72)  :         global_p_024_003                          
c     par( 73)  :         global_p_024_004                          
c     par( 74)  :         global_p_024_005                          
c     par( 75)  :         global_p_024_006                          
c     par( 76)  :         global_p_024_007                          
c     par( 77)  :         global_p_024_008                          
c     par( 78)  :         global_p_025_001                          
c     par( 79)  :         global_p_025_002                          
c     par( 80)  :         global_p_025_003                          
c     par( 81)  :         global_p_026_001                          
c     par( 82)  :         global_p_026_002                          
c     par( 83)  :         global_p_026_003                          
c     par( 84)  :         global_p_026_004                          
c     par( 85)  :         global_p_026_005                          
c     par( 86)  :         global_p_026_006                          
c     par( 87)  :         global_p_026_007                          
c     par( 88)  :         global_p_026_008                          
c     par( 89)  :         global_p_027_001                          
c     par( 90)  :         global_p_027_002                          
c     par( 91)  :         global_p_027_003                          
c     par( 92)  :         global_p_027_004                          
c     par( 93)  :         global_p_028_001                          
c     par( 94)  :         global_p_029_001                          
c     par( 95)  :         global_p_029_002                          
c     par( 96)  :         global_p_029_003                          
c     par( 97)  :         global_p_029_004                          
c     par( 98)  :         global_p_029_005                          
c     par( 99)  :         global_p_029_006                          
c     par(100)  :         global_p_030_001                          
c     par(101)  :         global_p_030_002                          
c     par(102)  :         global_p_030_003                          
c     par(103)  :         global_p_030_004                          
c     par(104)  :         global_p_030_005                          
c     par(105)  :         global_p_031_001                          
c     par(106)  :         global_p_031_002                          
c     par(107)  :         global_p_031_003                          
c     par(108)  :         global_p_032_001                          
c     par(109)  :         global_p_032_002                          
c     par(110)  :         global_p_033_001                          
c     par(111)  :         global_p_033_002                          
c     par(112)  :         global_p_033_003                          
c     par(113)  :         global_p_034_001                          
c     par(114)  :         global_p_034_002                          
c     par(115)  :         global_p_035_001                          
c     ---
c
      double precision    com(1)
      double precision    spe(33)
      double precision    par(115)
      double precision    rul(2)
      double precision    rea(57)
      logical             eve(0)
c
      common              /SBMLVARIABLES/ com, spe, par, rul, rea, eve
c
c
      double precision    fun1
      double precision    fun2
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
      rul(1) = par(94) * fun2(spe(25), par(95), par(96)) *
     &          (1.000000D+00 + par(97) * fun1(spe(24), par(98),
     &          par(99)))
      rul(2) = par(100) * (fun1(spe(24), par(101), par(102)) +
     &          fun2(spe(24), par(103), par(104)))
c 
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      rea(1) = (par(1) + par(2) * fun1(spe(24), par(3), par(4))) *
     &          fun2(spe(25), par(5), par(6))
      rea(2) = (par(7) + par(8) * fun1(spe(32), par(9), par(10))) *
     &          spe(1)
      rea(3) = (par(7) + par(8) * fun1(spe(32), par(9), par(10))) *
     &          spe(1) / par(11)
      rea(4) = par(12) * spe(2) * spe(3)
      rea(5) = par(13) * spe(2)
      rea(6) = par(14) * spe(5)
      rea(7) = par(15) * spe(4)
      rea(8) = par(16) / (1.000000D+00 + pow(spe(28) / par(17),
     &          par(18)) + pow(spe(27) / par(19), par(20))) *
     &          fun2(rul(1), par(21), par(22))
      rea(9) = (par(23) + par(24) * fun1(spe(32), par(25), par(26)))
     &          * spe(6)
      rea(10) = (par(23) + par(24) * fun1(spe(32), par(25),
     &          par(26))) * spe(6) / par(11)
      rea(11) = par(27) * spe(8) * spe(7)
      rea(12) = par(28) * spe(7)
      rea(13) = par(29) * spe(10)
      rea(14) = par(30) * spe(9)
      rea(15) = par(31) * fun1(spe(7), par(32), par(33))
      rea(16) = par(34) * fun1(spe(25), par(35), par(36)) * spe(11)
      rea(17) = par(38) * fun1(spe(9), par(39), par(40))
      rea(18) = par(41) * spe(9) * spe(12)
      rea(19) = par(42) * pow(spe(4) / par(43), par(44)) * spe(11) *
     &          spe(13)
      rea(20) = par(45) * spe(9) * spe(14) * (1.000000D+00 - spe(14)
     &          / par(46))
      rea(21) = par(47) * pow(spe(4) / par(43), par(48)) * spe(11) *
     &          spe(14)
      rea(22) = par(49) * pow(spe(4) / par(43), par(50)) * spe(15) *
     &          (1.000000D+00 - spe(15) / par(46))
      rea(23) = par(51) * (spe(4) / par(43)) * spe(11) * spe(15)
      rea(24) = par(52) * pow(spe(4) / par(43), par(53)) * spe(11) *
     &          spe(16)
      rea(25) = par(54) * pow(spe(4) / par(43), par(53)) * spe(11) *
     &          fun1(spe(16), par(55), par(56))
      rea(26) = par(57) * spe(17)
      rea(27) = par(58) * fun1(spe(17), par(59), par(60))
      rea(28) = par(61) * spe(18)
      rea(29) = par(62) * spe(19)
      rea(30) = par(63) * spe(20)
      rea(31) = par(63) * par(64) * fun1(spe(32), par(65), par(66))
     &          * spe(20)
      rea(32) = par(67) * spe(21)
      rea(33) = par(67) * par(64) * fun1(spe(32), par(65), par(66))
     &          * spe(21)
      rea(34) = par(68) * spe(22)
      rea(35) = par(68) * par(64) * fun1(spe(32), par(65), par(66))
     &          * spe(22)
      rea(36) = par(69) * (1.000000D+00 + par(64) * fun1(spe(32),
     &          par(65), par(66))) * spe(23)
      rea(37) = par(70) + par(71) * spe(13) + par(72) * spe(2) *
     &          spe(14) + par(73) * spe(15) + par(74) * spe(2) *
     &          spe(16) + par(75) * spe(20) + par(76) * spe(23)
      rea(38) = par(77) * spe(24)
      rea(39) = par(78) + par(79) * spe(23)
      rea(40) = par(80) * spe(25)
      rea(41) = par(81) + par(82) * spe(16) + par(83) * spe(18) +
     &          par(84) * spe(20) + par(85) * spe(21) + par(86) *
     &          spe(22) + par(87) * spe(23)
      rea(42) = par(88) * spe(26)
      rea(43) = par(89) + par(90) * spe(13) + par(91) * spe(19)
      rea(44) = par(92) * spe(27)
      rea(45) = par(93) * spe(28)
      rea(46) = rul(2) * rul(1)
      rea(47) = par(105) * spe(30) * spe(29)
      rea(48) = par(106) * spe(32)
      rea(49) = par(107) * spe(29)
      rea(50) = par(108) * spe(30)
      rea(51) = par(109) * spe(31)
      rea(52) = par(110)
      rea(53) = par(111) * spe(33)
      rea(54) = par(112) * spe(31)
      rea(55) = par(113) * spe(32)
      rea(56) = par(114) * spe(33)
      rea(57) = par(115) * spe(33)
c 
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      dy(1) = ( + rea(1) - rea(2) ) / com(1)
      dy(2) = ( + rea(3) - rea(4) - rea(5) ) / com(1)
      dy(3) = ( - rea(4) + rea(6) ) / com(1)
      dy(4) = ( + rea(4) - rea(7) ) / com(1)
      dy(5) = ( - rea(6) + rea(7) ) / com(1)
      dy(6) = ( + rea(8) - rea(9) ) / com(1)
      dy(7) = ( + rea(10) - rea(11) - rea(12) ) / com(1)
      dy(8) = ( - rea(11) + rea(13) ) / com(1)
      dy(9) = ( + rea(11) - rea(14) ) / com(1)
      dy(10) = ( - rea(13) + rea(14) ) / com(1)
      dy(11) = ( + rea(15) - rea(16) ) / com(1)
      dy(12) = ( + rea(17) - rea(18) ) / com(1)
      dy(13) = ( + rea(18) - rea(19) ) / com(1)
      dy(14) = ( + rea(19) + rea(20) - rea(21) ) / com(1)
      dy(15) = ( + rea(21) + rea(22) - rea(23) ) / com(1)
      dy(16) = ( + rea(23) - rea(24) ) / com(1)
      dy(17) = ( + rea(25) - rea(26) ) / com(1)
      dy(18) = ( + rea(27) - rea(28) ) / com(1)
      dy(19) = ( + rea(28) - rea(29) ) / com(1)
      dy(20) = ( + rea(29) - rea(30) - rea(31) ) / com(1)
      dy(21) = ( + rea(30) - rea(32) - rea(33) ) / com(1)
      dy(22) = ( + rea(32) - rea(34) - rea(35) ) / com(1)
      dy(23) = ( + rea(34) - rea(36) ) / com(1)
      dy(24) = ( + rea(37) - rea(38) ) / com(1)
      dy(25) = ( + rea(39) - rea(40) ) / com(1)
      dy(26) = ( + rea(41) - rea(42) ) / com(1)
      dy(27) = ( + rea(43) - rea(44) ) / com(1)
      dy(28) = ( + rea(42) - rea(45) ) / com(1)
      dy(29) = ( + rea(46) - rea(47) + rea(48) - rea(49) ) / com(1)
      dy(30) = ( - rea(47) + rea(48) - rea(50) + rea(51) ) / com(1)
      dy(31) = ( + rea(50) - rea(51) + rea(52) + rea(53) - rea(54) )
     &          / com(1)
      dy(32) = ( + rea(47) - rea(48) - rea(55) + rea(56) ) / com(1)
      dy(33) = ( - rea(53) + rea(55) - rea(56) - rea(57) ) / com(1)
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
      double precision    com(1)
      double precision    spe(33)
      double precision    par(115)
      double precision    rul(2)
      double precision    rea(57)
      logical             eve(0)
c
      common              /SBMLVARIABLES/ com, spe, par, rul, rea, eve
c
c
      double precision    fun1
      double precision    fun2
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
c     ---
c
c     ---
c     spe(  1)  :         LH_pit                                    
c     spe(  2)  :         LH_blood                                  
c     spe(  3)  :         R_LH                                      
c     spe(  4)  :         LH_R                                      
c     spe(  5)  :         R_LH_des                                  
c     spe(  6)  :         FSH_pit                                   
c     spe(  7)  :         FSH_blood                                 
c     spe(  8)  :         R_FSH                                     
c     spe(  9)  :         FSH_R                                     
c     spe( 10)  :         R_FSH_des                                 
c     spe( 11)  :         S_foll                                    
c     spe( 12)  :         AF1                                       
c     spe( 13)  :         AF2                                       
c     spe( 14)  :         AF3                                       
c     spe( 15)  :         AF4                                       
c     spe( 16)  :         PrF                                       
c     spe( 17)  :         OvF                                       
c     spe( 18)  :         Sc1                                       
c     spe( 19)  :         Sc2                                       
c     spe( 20)  :         Lut1                                      
c     spe( 21)  :         Lut2                                      
c     spe( 22)  :         Lut3                                      
c     spe( 23)  :         Lut4                                      
c     spe( 24)  :         E2                                        
c     spe( 25)  :         P4                                        
c     spe( 26)  :         InhA                                      
c     spe( 27)  :         InhB                                      
c     spe( 28)  :         InhA_tau                                  
c     spe( 29)  :         GnRH_G                                    
c     spe( 30)  :         R_GnRH_a                                  
c     spe( 31)  :         R_GnRH_i                                  
c     spe( 32)  :         GnRH_R_a                                  
c     spe( 33)  :         GnRH_R_i                                  
c     ---
c
c     ---
c     par(  1)  :         global_p_001_001                          
c     par(  2)  :         global_p_001_002                          
c     par(  3)  :         global_p_001_003                          
c     par(  4)  :         global_p_001_004                          
c     par(  5)  :         global_p_001_005                          
c     par(  6)  :         global_p_001_006                          
c     par(  7)  :         global_p_001_007                          
c     par(  8)  :         global_p_001_008                          
c     par(  9)  :         global_p_001_009                          
c     par( 10)  :         global_p_001_010                          
c     par( 11)  :         global_p_002_001                          
c     par( 12)  :         global_p_002_002                          
c     par( 13)  :         global_p_002_003                          
c     par( 14)  :         global_p_003_001                          
c     par( 15)  :         global_p_004_001                          
c     par( 16)  :         global_p_006_001                          
c     par( 17)  :         global_p_006_002                          
c     par( 18)  :         global_p_006_003                          
c     par( 19)  :         global_p_006_004                          
c     par( 20)  :         global_p_006_005                          
c     par( 21)  :         global_p_006_006                          
c     par( 22)  :         global_p_006_007                          
c     par( 23)  :         global_p_006_008                          
c     par( 24)  :         global_p_006_009                          
c     par( 25)  :         global_p_006_010                          
c     par( 26)  :         global_p_006_011                          
c     par( 27)  :         global_p_007_001                          
c     par( 28)  :         global_p_007_002                          
c     par( 29)  :         global_p_008_001                          
c     par( 30)  :         global_p_009_001                          
c     par( 31)  :         global_p_011_001                          
c     par( 32)  :         global_p_011_002                          
c     par( 33)  :         global_p_011_003                          
c     par( 34)  :         global_p_011_004                          
c     par( 35)  :         global_p_011_005                          
c     par( 36)  :         global_p_011_006                          
c     par( 37)  :         global_p_011_007                          
c     par( 38)  :         global_p_012_001                          
c     par( 39)  :         global_p_012_002                          
c     par( 40)  :         global_p_012_003                          
c     par( 41)  :         global_p_012_004                          
c     par( 42)  :         global_p_013_001                          
c     par( 43)  :         global_p_013_002                          
c     par( 44)  :         global_p_013_003                          
c     par( 45)  :         global_p_014_001                          
c     par( 46)  :         global_p_014_002                          
c     par( 47)  :         global_p_014_003                          
c     par( 48)  :         global_p_014_004                          
c     par( 49)  :         global_p_015_001                          
c     par( 50)  :         global_p_015_002                          
c     par( 51)  :         global_p_015_003                          
c     par( 52)  :         global_p_016_001                          
c     par( 53)  :         global_p_016_002                          
c     par( 54)  :         global_p_017_001                          
c     par( 55)  :         global_p_017_002                          
c     par( 56)  :         global_p_017_003                          
c     par( 57)  :         global_p_017_004                          
c     par( 58)  :         global_p_018_001                          
c     par( 59)  :         global_p_018_002                          
c     par( 60)  :         global_p_018_003                          
c     par( 61)  :         global_p_018_004                          
c     par( 62)  :         global_p_019_001                          
c     par( 63)  :         global_p_020_001                          
c     par( 64)  :         global_p_020_002                          
c     par( 65)  :         global_p_020_003                          
c     par( 66)  :         global_p_020_004                          
c     par( 67)  :         global_p_021_001                          
c     par( 68)  :         global_p_022_001                          
c     par( 69)  :         global_p_023_001                          
c     par( 70)  :         global_p_024_001                          
c     par( 71)  :         global_p_024_002                          
c     par( 72)  :         global_p_024_003                          
c     par( 73)  :         global_p_024_004                          
c     par( 74)  :         global_p_024_005                          
c     par( 75)  :         global_p_024_006                          
c     par( 76)  :         global_p_024_007                          
c     par( 77)  :         global_p_024_008                          
c     par( 78)  :         global_p_025_001                          
c     par( 79)  :         global_p_025_002                          
c     par( 80)  :         global_p_025_003                          
c     par( 81)  :         global_p_026_001                          
c     par( 82)  :         global_p_026_002                          
c     par( 83)  :         global_p_026_003                          
c     par( 84)  :         global_p_026_004                          
c     par( 85)  :         global_p_026_005                          
c     par( 86)  :         global_p_026_006                          
c     par( 87)  :         global_p_026_007                          
c     par( 88)  :         global_p_026_008                          
c     par( 89)  :         global_p_027_001                          
c     par( 90)  :         global_p_027_002                          
c     par( 91)  :         global_p_027_003                          
c     par( 92)  :         global_p_027_004                          
c     par( 93)  :         global_p_028_001                          
c     par( 94)  :         global_p_029_001                          
c     par( 95)  :         global_p_029_002                          
c     par( 96)  :         global_p_029_003                          
c     par( 97)  :         global_p_029_004                          
c     par( 98)  :         global_p_029_005                          
c     par( 99)  :         global_p_029_006                          
c     par(100)  :         global_p_030_001                          
c     par(101)  :         global_p_030_002                          
c     par(102)  :         global_p_030_003                          
c     par(103)  :         global_p_030_004                          
c     par(104)  :         global_p_030_005                          
c     par(105)  :         global_p_031_001                          
c     par(106)  :         global_p_031_002                          
c     par(107)  :         global_p_031_003                          
c     par(108)  :         global_p_032_001                          
c     par(109)  :         global_p_032_002                          
c     par(110)  :         global_p_033_001                          
c     par(111)  :         global_p_033_002                          
c     par(112)  :         global_p_033_003                          
c     par(113)  :         global_p_034_001                          
c     par(114)  :         global_p_034_002                          
c     par(115)  :         global_p_035_001                          
c     ---
c
c-----------------------------------------------------------------------
c
      double precision    com(1)
      double precision    spe(33)
      double precision    par(115)
      double precision    rul(2)
      double precision    rea(57)
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
c 
c
      spe(  1) = 2.596960D+05
      spe(  2) = 2.710000D+00
      spe(  3) = 8.398000D+00
      spe(  4) = 2.660000D-01
      spe(  5) = 7.100000D-01
      spe(  6) = 5.973400D+04
      spe(  7) = 5.280000D+00
      spe(  8) = 5.903000D+00
      spe(  9) = 7.960000D-01
      spe( 10) = 1.800000D+00
      spe( 11) = 1.560000D-01
      spe( 12) = 2.593000D+00
      spe( 13) = 2.314000D+01
      spe( 14) = 4.920000D-01
      spe( 15) = 1.610000D-05
      spe( 16) = 2.449000D-01
      spe( 17) = 0.000000D+00
      spe( 18) = 1.000000D-08
      spe( 19) = 2.000000D-06
      spe( 20) = 2.000000D-05
      spe( 21) = 3.000000D-04
      spe( 22) = 3.000000D-03
      spe( 23) = 1.100000D-02
      spe( 24) = 4.136000D+01
      spe( 25) = 1.978000D+00
      spe( 26) = 9.300000D-01
      spe( 27) = 6.053000D+01
      spe( 28) = 7.114000D+01
      spe( 29) = 1.520000D-02
      spe( 30) = 9.000000D-03
      spe( 31) = 9.600000D-04
      spe( 32) = 6.500000D-05
      spe( 33) = 5.900000D-05
c 
c
      par(  1) = 7.309920D+03
      par(  2) = 7.309920D+03
      par(  3) = 1.922000D+02
      par(  4) = 1.000000D+01
      par(  5) = 2.371000D+00
      par(  6) = 1.000000D+00
      par(  7) = 4.760000D-03
      par(  8) = 1.904000D-01
      par(  9) = 3.000000D-04
      par( 10) = 5.000000D+00
      par( 11) = 5.000000D+00
      par( 12) = 2.143000D+00
      par( 13) = 7.485100D+01
      par( 14) = 6.894900D+01
      par( 15) = 1.833600D+02
      par( 16) = 2.213000D+04
      par( 17) = 9.581000D+01
      par( 18) = 5.000000D+00
      par( 19) = 7.000000D+01
      par( 20) = 3.000000D+00
      par( 21) = 1.000000D+01
      par( 22) = 3.000000D+00
      par( 23) = 5.700000D-02
      par( 24) = 2.720000D-01
      par( 25) = 3.000000D-04
      par( 26) = 3.000000D+00
      par( 27) = 3.529000D+00
      par( 28) = 1.142500D+02
      par( 29) = 6.102900D+01
      par( 30) = 1.383000D+02
      par( 31) = 2.190000D-01
      par( 32) = 3.000000D+00
      par( 33) = 5.000000D+00
      par( 34) = 1.343000D+00
      par( 35) = 1.235000D+00
      par( 36) = 5.000000D+00
      par( 37) = 2.000000D-01
      par( 38) = 3.663000D+00
      par( 39) = 6.080000D-01
      par( 40) = 3.000000D+00
      par( 41) = 1.221000D+00
      par( 42) = 4.882000D+00
      par( 43) = 2.726000D+00
      par( 44) = 3.689000D+00
      par( 45) = 1.220000D-01
      par( 46) = 1.000000D+01
      par( 47) = 1.220600D+02
      par( 48) = 5.000000D+00
      par( 49) = 1.220600D+01
      par( 50) = 2.000000D+00
      par( 51) = 3.327500D+02
      par( 52) = 1.220600D+02
      par( 53) = 6.000000D+00
      par( 54) = 7.984000D+00
      par( 55) = 3.000000D+00
      par( 56) = 1.000000D+01
      par( 57) = 1.220600D+01
      par( 58) = 1.208000D+00
      par( 59) = 2.000000D-02
      par( 60) = 1.000000D+01
      par( 61) = 1.221000D+00
      par( 62) = 9.580000D-01
      par( 63) = 9.250000D-01
      par( 64) = 2.000000D+01
      par( 65) = 8.000000D-04
      par( 66) = 5.000000D+00
      par( 67) = 7.567000D-01
      par( 68) = 6.100000D-01
      par( 69) = 5.430000D-01
      par( 70) = 5.155800D+01
      par( 71) = 2.094500D+00
      par( 72) = 9.280000D+00
      par( 73) = 3.480270D+03
      par( 74) = 9.720000D-01
      par( 75) = 1.713710D+03
      par( 76) = 8.675140D+03
      par( 77) = 5.235000D+00
      par( 78) = 9.430000D-01
      par( 79) = 7.616400D+02
      par( 80) = 5.130000D+00
      par( 81) = 1.445000D+00
      par( 82) = 2.285000D+00
      par( 83) = 6.000000D+01
      par( 84) = 1.800000D+02
      par( 85) = 2.821100D+01
      par( 86) = 1.940700D+02
      par( 87) = 1.142500D+02
      par( 88) = 4.287000D+00
      par( 89) = 8.994300D+01
      par( 90) = 4.474700D+02
      par( 91) = 1.322402D+05
      par( 92) = 1.724500D+02
      par( 93) = 1.990000D-01
      par( 94) = 1.600000D+01
      par( 95) = 3.000000D+00
      par( 96) = 2.000000D+00
      par( 97) = 1.000000D+00
      par( 98) = 2.200000D+02
      par( 99) = 1.000000D+01
      par(100) = 5.593000D-03
      par(101) = 2.200000D+02
      par(102) = 2.000000D+00
      par(103) = 9.600000D+00
      par(104) = 1.000000D+00
      par(105) = 3.221800D+02
      par(106) = 6.443500D+02
      par(107) = 4.470000D-01
      par(108) = 3.222000D+00
      par(109) = 3.221800D+01
      par(110) = 8.949000D-05
      par(111) = 3.221800D+01
      par(112) = 8.950000D-02
      par(113) = 3.221800D+01
      par(114) = 3.222000D+00
      par(115) = 8.950000D-03
c 
c
c 
c
c-----------------------------------------------------------------------
c
      mc = min0(1,lc)
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
         lc = 1
c
      end if
c
c-----------------------------------------------------------------------
c
      ms = min0(33,ls)
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
         ls = 33
c
      end if
c
c-----------------------------------------------------------------------
c
      mp = min0(115,lp)
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
         lp = 115
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
      character           idc(1)*(*)
c
      integer             ncom
c
c-----------------------------------------------------------------------
c
      ncom = 1
c     ---
      idc(  1) = 'default'//char(0)
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
      character           ids(33)*(*)
c
      integer             j, nspe
c
c-----------------------------------------------------------------------
c
      nspe = 33
c     ---
      do j = 1, nspe
c
        ids(j) = char(0)
c
      end do
c     ---
      ids(  1) = 'LH_pit'//char(0)
      ids(  2) = 'LH_blood'//char(0)
      ids(  3) = 'R_LH'//char(0)
      ids(  4) = 'LH_R'//char(0)
      ids(  5) = 'R_LH_des'//char(0)
      ids(  6) = 'FSH_pit'//char(0)
      ids(  7) = 'FSH_blood'//char(0)
      ids(  8) = 'R_FSH'//char(0)
      ids(  9) = 'FSH_R'//char(0)
      ids( 10) = 'R_FSH_des'//char(0)
      ids( 11) = 'S_foll'//char(0)
      ids( 12) = 'AF1'//char(0)
      ids( 13) = 'AF2'//char(0)
      ids( 14) = 'AF3'//char(0)
      ids( 15) = 'AF4'//char(0)
      ids( 16) = 'PrF'//char(0)
      ids( 17) = 'OvF'//char(0)
      ids( 18) = 'Sc1'//char(0)
      ids( 19) = 'Sc2'//char(0)
      ids( 20) = 'Lut1'//char(0)
      ids( 21) = 'Lut2'//char(0)
      ids( 22) = 'Lut3'//char(0)
      ids( 23) = 'Lut4'//char(0)
      ids( 24) = 'E2'//char(0)
      ids( 25) = 'P4'//char(0)
      ids( 26) = 'InhA'//char(0)
      ids( 27) = 'InhB'//char(0)
      ids( 28) = 'InhA_tau'//char(0)
      ids( 29) = 'GnRH_G'//char(0)
      ids( 30) = 'R_GnRH_a'//char(0)
      ids( 31) = 'R_GnRH_i'//char(0)
      ids( 32) = 'GnRH_R_a'//char(0)
      ids( 33) = 'GnRH_R_i'//char(0)
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
      character           idp(115)*(*)
c
      integer             j, npar
c
c-----------------------------------------------------------------------
c
      npar = 115
c     ---
      do j = 1, npar
c
        idp(j) = char(0)
c
      end do
c     ---
      idp(  1) = 'global_p_001_001'//char(0)
      idp(  2) = 'global_p_001_002'//char(0)
      idp(  3) = 'global_p_001_003'//char(0)
      idp(  4) = 'global_p_001_004'//char(0)
      idp(  5) = 'global_p_001_005'//char(0)
      idp(  6) = 'global_p_001_006'//char(0)
      idp(  7) = 'global_p_001_007'//char(0)
      idp(  8) = 'global_p_001_008'//char(0)
      idp(  9) = 'global_p_001_009'//char(0)
      idp( 10) = 'global_p_001_010'//char(0)
      idp( 11) = 'global_p_002_001'//char(0)
      idp( 12) = 'global_p_002_002'//char(0)
      idp( 13) = 'global_p_002_003'//char(0)
      idp( 14) = 'global_p_003_001'//char(0)
      idp( 15) = 'global_p_004_001'//char(0)
      idp( 16) = 'global_p_006_001'//char(0)
      idp( 17) = 'global_p_006_002'//char(0)
      idp( 18) = 'global_p_006_003'//char(0)
      idp( 19) = 'global_p_006_004'//char(0)
      idp( 20) = 'global_p_006_005'//char(0)
      idp( 21) = 'global_p_006_006'//char(0)
      idp( 22) = 'global_p_006_007'//char(0)
      idp( 23) = 'global_p_006_008'//char(0)
      idp( 24) = 'global_p_006_009'//char(0)
      idp( 25) = 'global_p_006_010'//char(0)
      idp( 26) = 'global_p_006_011'//char(0)
      idp( 27) = 'global_p_007_001'//char(0)
      idp( 28) = 'global_p_007_002'//char(0)
      idp( 29) = 'global_p_008_001'//char(0)
      idp( 30) = 'global_p_009_001'//char(0)
      idp( 31) = 'global_p_011_001'//char(0)
      idp( 32) = 'global_p_011_002'//char(0)
      idp( 33) = 'global_p_011_003'//char(0)
      idp( 34) = 'global_p_011_004'//char(0)
      idp( 35) = 'global_p_011_005'//char(0)
      idp( 36) = 'global_p_011_006'//char(0)
      idp( 37) = 'global_p_011_007'//char(0)
      idp( 38) = 'global_p_012_001'//char(0)
      idp( 39) = 'global_p_012_002'//char(0)
      idp( 40) = 'global_p_012_003'//char(0)
      idp( 41) = 'global_p_012_004'//char(0)
      idp( 42) = 'global_p_013_001'//char(0)
      idp( 43) = 'global_p_013_002'//char(0)
      idp( 44) = 'global_p_013_003'//char(0)
      idp( 45) = 'global_p_014_001'//char(0)
      idp( 46) = 'global_p_014_002'//char(0)
      idp( 47) = 'global_p_014_003'//char(0)
      idp( 48) = 'global_p_014_004'//char(0)
      idp( 49) = 'global_p_015_001'//char(0)
      idp( 50) = 'global_p_015_002'//char(0)
      idp( 51) = 'global_p_015_003'//char(0)
      idp( 52) = 'global_p_016_001'//char(0)
      idp( 53) = 'global_p_016_002'//char(0)
      idp( 54) = 'global_p_017_001'//char(0)
      idp( 55) = 'global_p_017_002'//char(0)
      idp( 56) = 'global_p_017_003'//char(0)
      idp( 57) = 'global_p_017_004'//char(0)
      idp( 58) = 'global_p_018_001'//char(0)
      idp( 59) = 'global_p_018_002'//char(0)
      idp( 60) = 'global_p_018_003'//char(0)
      idp( 61) = 'global_p_018_004'//char(0)
      idp( 62) = 'global_p_019_001'//char(0)
      idp( 63) = 'global_p_020_001'//char(0)
      idp( 64) = 'global_p_020_002'//char(0)
      idp( 65) = 'global_p_020_003'//char(0)
      idp( 66) = 'global_p_020_004'//char(0)
      idp( 67) = 'global_p_021_001'//char(0)
      idp( 68) = 'global_p_022_001'//char(0)
      idp( 69) = 'global_p_023_001'//char(0)
      idp( 70) = 'global_p_024_001'//char(0)
      idp( 71) = 'global_p_024_002'//char(0)
      idp( 72) = 'global_p_024_003'//char(0)
      idp( 73) = 'global_p_024_004'//char(0)
      idp( 74) = 'global_p_024_005'//char(0)
      idp( 75) = 'global_p_024_006'//char(0)
      idp( 76) = 'global_p_024_007'//char(0)
      idp( 77) = 'global_p_024_008'//char(0)
      idp( 78) = 'global_p_025_001'//char(0)
      idp( 79) = 'global_p_025_002'//char(0)
      idp( 80) = 'global_p_025_003'//char(0)
      idp( 81) = 'global_p_026_001'//char(0)
      idp( 82) = 'global_p_026_002'//char(0)
      idp( 83) = 'global_p_026_003'//char(0)
      idp( 84) = 'global_p_026_004'//char(0)
      idp( 85) = 'global_p_026_005'//char(0)
      idp( 86) = 'global_p_026_006'//char(0)
      idp( 87) = 'global_p_026_007'//char(0)
      idp( 88) = 'global_p_026_008'//char(0)
      idp( 89) = 'global_p_027_001'//char(0)
      idp( 90) = 'global_p_027_002'//char(0)
      idp( 91) = 'global_p_027_003'//char(0)
      idp( 92) = 'global_p_027_004'//char(0)
      idp( 93) = 'global_p_028_001'//char(0)
      idp( 94) = 'global_p_029_001'//char(0)
      idp( 95) = 'global_p_029_002'//char(0)
      idp( 96) = 'global_p_029_003'//char(0)
      idp( 97) = 'global_p_029_004'//char(0)
      idp( 98) = 'global_p_029_005'//char(0)
      idp( 99) = 'global_p_029_006'//char(0)
      idp(100) = 'global_p_030_001'//char(0)
      idp(101) = 'global_p_030_002'//char(0)
      idp(102) = 'global_p_030_003'//char(0)
      idp(103) = 'global_p_030_004'//char(0)
      idp(104) = 'global_p_030_005'//char(0)
      idp(105) = 'global_p_031_001'//char(0)
      idp(106) = 'global_p_031_002'//char(0)
      idp(107) = 'global_p_031_003'//char(0)
      idp(108) = 'global_p_032_001'//char(0)
      idp(109) = 'global_p_032_002'//char(0)
      idp(110) = 'global_p_033_001'//char(0)
      idp(111) = 'global_p_033_002'//char(0)
      idp(112) = 'global_p_033_003'//char(0)
      idp(113) = 'global_p_034_001'//char(0)
      idp(114) = 'global_p_034_002'//char(0)
      idp(115) = 'global_p_035_001'//char(0)
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
      idm(  1) = 'PAEON_V2'//char(0)
      idm(  2) = 'Mon Feb  2 18:21:20 2015'//char(0)
      idm(  3) = '001422897680'//char(0)
      idm(  4) = 'PAEON_V2.xml'//char(0)
c     ---
c
c-----------------------------------------------------------------------
c
      return
      end subroutine get_model_ids
c
c=======================================================================
c
c=======================================================================
c
      double precision function fun1 ( x1, x2, x3 )
c 
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x1
      double precision    x2
      double precision    x3
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
c
      double precision    pi
      double precision    eul
c
      parameter           ( pi  = 3.141592653589793238462643383276D0 )
      parameter           ( eul = 2.718281828459045235360287471352D0 )
c
c-----------------------------------------------------------------------
c
      fun1 = pow(x1 / x2, x3) / (1.000000D+00 + pow(x1 / x2, x3))
c 
c-----------------------------------------------------------------------
c
      return
      end function fun1
c
c=======================================================================
c
      double precision function fun2 ( x1, x2, x3 )
c 
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x1
      double precision    x2
      double precision    x3
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
c
      double precision    pi
      double precision    eul
c
      parameter           ( pi  = 3.141592653589793238462643383276D0 )
      parameter           ( eul = 2.718281828459045235360287471352D0 )
c
c-----------------------------------------------------------------------
c
      fun2 = 1.000000D+00 / (1.000000D+00 + pow(x1 / x2, x3))
c 
c-----------------------------------------------------------------------
c
      return
      end function fun2
c
c 
c
c 
c

