c-----
c SBML Model : Hettling2011                                           
c              ~~~~~~~~~~~~
c       Date : Thu Sep 11 10:11:31 2014
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
c     com(  1)  :         IMS                                       
c     com(  2)  :         CYT                                       
c     com(  3)  :         cell                                      
c     ---
c
c     ---
c     spe(  1)  :         ADPi                                      
c     spe(  2)  :         ATPi                                      
c     spe(  3)  :         Cri                                       
c     spe(  4)  :         PCri                                      
c     spe(  5)  :         PCr                                       
c     spe(  6)  :         ADP                                       
c     spe(  7)  :         ATP                                       
c     spe(  8)  :         Cr                                        
c     spe(  9)  :         P_ii                                      
c     spe( 10)  :         P_i                                       
c     ---
c
c     ---
c     par(  1)  :         global_j_diff_pcr                         
c     par(  2)  :         global_j_diff_atp                         
c     par(  3)  :         global_densyn                             
c     par(  4)  :         global_tmito                              
c     par(  5)  :         global_vatpnorm                           
c     par(  6)  :         global_jsyn                               
c     par(  7)  :         global_r_diff_pcr                         
c     par(  8)  :         global_j_ck_mi                            
c     par(  9)  :         global_j_ck_mm                            
c     par( 10)  :         global_j_diff_adp                         
c     par( 11)  :         global_j_diff_cr                          
c     par( 12)  :         global_j_diff_pi                          
c     par( 13)  :         global_stepsize                           
c     par( 14)  :         global_phase                              
c     par( 15)  :         global_heartrate_bpm                      
c     par( 16)  :         global_heartrate_basis                    
c     par( 17)  :         global_heartrate_test                     
c     par( 18)  :         global_fracDia                            
c     par( 19)  :         global_fracSysUp                          
c     par( 20)  :         global_fracSysDown                        
c     par( 21)  :         global_VhydAmp_basis                      
c     par( 22)  :         global_VhydAmp_test                       
c     par( 23)  :         global_time_Jhyd_step                     
c     par( 24)  :         global_Jhyd_test                          
c     par( 25)  :         global_Jhyd_basis                         
c     par( 26)  :         global_last_time_fired                    
c     par( 27)  :         global_Jhyd                               
c     par( 28)  :         global_ck_factor_iaa                      
c     par( 29)  :         global_ck_factor_ia                       
c     par( 30)  :         global_tmito_factor                       
c     par( 31)  :         global_PSmomATP                           
c     par( 32)  :         global_K_CK_eq                            
c     par( 33)  :         global_VmaxMMb                            
c     par( 34)  :         global_VmaxMib                            
c     par( 35)  :         global_VmaxMif_full_activity              
c     par( 36)  :         global_VmaxMMf_full_activity              
c     par( 37)  :         global_VmaxMif                            
c     par( 38)  :         global_VmaxMMf                            
c     par( 39)  :         global_KiaMi                              
c     par( 40)  :         global_KbMi                               
c     par( 41)  :         global_KicMi                              
c     par( 42)  :         global_KdMi                               
c     par( 43)  :         global_KibMi                              
c     par( 44)  :         global_KidMi                              
c     par( 45)  :         global_KiaMM                              
c     par( 46)  :         global_KbMM                               
c     par( 47)  :         global_KicMM                              
c     par( 48)  :         global_KdMM                               
c     par( 49)  :         global_KibMM                              
c     par( 50)  :         global_KidMM                              
c     par( 51)  :         global_Vmaxsyn                            
c     par( 52)  :         global_Kadp                               
c     par( 53)  :         global_Kpi                                
c     par( 54)  :         global_PSmomPi                            
c     par( 55)  :         global_PSmomCr                            
c     par( 56)  :         global_PSmomPCr                           
c     par( 57)  :         global_pulsatility                        
c     ---
c
      double precision    com(3)
      double precision    spe(10)
      double precision    par(57)
      double precision    rul(20)
      double precision    rea(9)
      logical             eve(3)
c
      common              /SBMLVARIABLES/ com, spe, par, rul, rea, eve
c
c
c 
      double precision    piecewise2
      double precision    piecewise6
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
c
      rul(1) = par(5)
      rul(2) = 1.000000D+00 + spe(1) / par(52) + spe(9) / par(53) +
     &          spe(1) * spe(9) / (par(52) * par(53))
      rul(3) = par(51) * (spe(1) * spe(9) / (par(53) * par(52) *
     &          rul(2)))
      rul(4) = par(30) * ((par(24) - rul(3)) / (par(24) - par(25)))
      rul(5) = par(35) * par(28)
      rul(6) = par(36) * par(28)
      rul(7) = par(32) * par(41) * par(42) * rul(5) / (par(39) *
     &          par(40))
      rul(8) = par(32) * par(47) * par(48) * rul(6) / (par(45) *
     &          par(46))
      rul(9) = par(56) * (spe(4) - spe(5))
      rul(10) = par(31) * (spe(2) - spe(7))
      rul(11) = rul(9) / (rul(9) + rul(10))
      rul(12) = (abs(t - par(23)) - floor(abs(t - par(23)) /
     &          (6.000000D+01 / par(15))) * 6.000000D+01 / par(15)) /
     &          (6.000000D+01 / par(15))
      rul(13) = 1.000000D+00 - par(19) - par(20)
      rul(14) = 2.000000D+00 * par(25) / (par(19) + par(20))
      rul(15) = 2.000000D+00 * par(24) / (par(19) + par(20))
      rul(16) = (rul(5) * spe(2) * spe(3) / (par(39) * par(40)) -
     &          rul(7) * spe(1) * spe(4) / (par(41) * par(42))) /
     &          (1.000000D+00 + spe(3) / par(43) + spe(4) / par(44) +
     &          spe(2) * (1.000000D+00 / par(39) + spe(3) / (par(39)
     &          * par(40))) + spe(1) * (1.000000D+00 / par(41) +
     &          spe(3) / (par(41) * par(43)) + spe(4) / (par(44) *
     &          (par(41) * par(42) / par(44)))))
      rul(17) = (rul(6) * spe(7) * spe(8) / (par(45) * par(46)) -
     &          rul(8) * spe(6) * spe(5) / (par(47) * par(48))) /
     &          (1.000000D+00 + spe(8) / par(49) + spe(5) / par(50) +
     &          spe(7) * (1.000000D+00 / par(45) + spe(8) / (par(45)
     &          * par(46))) + spe(6) * (1.000000D+00 / par(47) +
     &          spe(8) / (par(47) * par(49)) + spe(5) / (par(50) *
     &          (par(47) * par(48) / par(50)))))
      rul(18) = par(31) * (spe(1) - spe(6))
      rul(19) = par(54) * (spe(9) - spe(10))
      rul(20) = par(55) * (spe(3) - spe(8))
c 
c
      rea(1) = par(51) * spe(1) * spe(9) / (par(52) * par(53) *
     &          (1.000000D+00 + spe(1) / par(52) + spe(9) / par(53) +
     &          spe(1) * spe(9) / (par(52) * par(53))))
      rea(2) = rul(16)
      rea(3) = rul(17)
      rea(4) = par(27)
      rea(5) = rul(19)
      rea(6) = rul(20)
      rea(7) = rul(18)
      rea(8) = rul(9)
      rea(9) = rul(10)
c 
c
c-----------------------------------------------------------------------
c
c
      eve(1) = and(and(geq(t, par(23)), gt(t - par(26), par(13))),
     &          eq(par(57), 1.000000D+00))
c 
      if ( eve(1) ) then
c
      par(27) = piecewise6((1.000000D+00 - (rul(12) - par(19)) /
     &          par(20)) * rul(15), and(gt(rul(12), par(19)),
     &          leq(rul(12), 1.000000D+00 - rul(13))), rul(12) /
     &          par(19) * rul(15), leq(rul(12), par(19)),
     &          0.000000D+00, geq(rul(12), 1.000000D+00 - rul(13)),
     &          par(27))
      par(26) = t
c 
      end if
c
c
      eve(2) = and(and(lt(t, par(23)), geq(t - par(26), par(13))),
     &          eq(par(57), 1.000000D+00))
c 
      if ( eve(2) ) then
c
      par(27) = piecewise6(0.000000D+00, leq(rul(12), rul(13)),
     &          (rul(12) - rul(13)) / par(20) * rul(14),
     &          and(gt(rul(12), rul(13)), leq(rul(12), 1.000000D+00 -
     &          par(19))), (1.000000D+00 - rul(12)) * rul(14) /
     &          par(19), gt(rul(12), 1.000000D+00 - par(19)),
     &          par(27))
      par(26) = t
c 
      end if
c
c
      eve(3) = geq(t, par(23))
c 
      if ( eve(3) ) then
c
      par(27) = piecewise2(par(24), eq(par(57), 0.000000D+00),
     &          par(27))
      par(15) = par(17)
      par(30) = 1.000000D+00
c 
      end if
c
c 
c-----------------------------------------------------------------------
c
      dy(1) = ( - rea(1) + rea(2) - rea(7) ) / com(1)
      dy(2) = ( + rea(1) - rea(2) - rea(9) ) / com(1)
      dy(3) = ( - rea(2) - rea(6) ) / com(1)
      dy(4) = ( + rea(2) - rea(8) ) / com(1)
      dy(5) = ( + rea(3) + rea(8) ) / com(2)
      dy(6) = ( + rea(3) + rea(4) + rea(7) ) / com(2)
      dy(7) = ( - rea(3) - rea(4) + rea(9) ) / com(2)
      dy(8) = ( - rea(3) + rea(6) ) / com(2)
      dy(9) = ( - rea(1) - rea(5) ) / com(1)
      dy(10) = ( + rea(4) + rea(5) ) / com(2)
c 
c
c-----------------------------------------------------------------------
c
      return
      end subroutine ydot
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
c     com(  1)  :         IMS                                       
c     com(  2)  :         CYT                                       
c     com(  3)  :         cell                                      
c     ---
c
c     ---
c     spe(  1)  :         ADPi                                      
c     spe(  2)  :         ATPi                                      
c     spe(  3)  :         Cri                                       
c     spe(  4)  :         PCri                                      
c     spe(  5)  :         PCr                                       
c     spe(  6)  :         ADP                                       
c     spe(  7)  :         ATP                                       
c     spe(  8)  :         Cr                                        
c     spe(  9)  :         P_ii                                      
c     spe( 10)  :         P_i                                       
c     ---
c
c     ---
c     par(  1)  :         global_j_diff_pcr                         
c     par(  2)  :         global_j_diff_atp                         
c     par(  3)  :         global_densyn                             
c     par(  4)  :         global_tmito                              
c     par(  5)  :         global_vatpnorm                           
c     par(  6)  :         global_jsyn                               
c     par(  7)  :         global_r_diff_pcr                         
c     par(  8)  :         global_j_ck_mi                            
c     par(  9)  :         global_j_ck_mm                            
c     par( 10)  :         global_j_diff_adp                         
c     par( 11)  :         global_j_diff_cr                          
c     par( 12)  :         global_j_diff_pi                          
c     par( 13)  :         global_stepsize                           
c     par( 14)  :         global_phase                              
c     par( 15)  :         global_heartrate_bpm                      
c     par( 16)  :         global_heartrate_basis                    
c     par( 17)  :         global_heartrate_test                     
c     par( 18)  :         global_fracDia                            
c     par( 19)  :         global_fracSysUp                          
c     par( 20)  :         global_fracSysDown                        
c     par( 21)  :         global_VhydAmp_basis                      
c     par( 22)  :         global_VhydAmp_test                       
c     par( 23)  :         global_time_Jhyd_step                     
c     par( 24)  :         global_Jhyd_test                          
c     par( 25)  :         global_Jhyd_basis                         
c     par( 26)  :         global_last_time_fired                    
c     par( 27)  :         global_Jhyd                               
c     par( 28)  :         global_ck_factor_iaa                      
c     par( 29)  :         global_ck_factor_ia                       
c     par( 30)  :         global_tmito_factor                       
c     par( 31)  :         global_PSmomATP                           
c     par( 32)  :         global_K_CK_eq                            
c     par( 33)  :         global_VmaxMMb                            
c     par( 34)  :         global_VmaxMib                            
c     par( 35)  :         global_VmaxMif_full_activity              
c     par( 36)  :         global_VmaxMMf_full_activity              
c     par( 37)  :         global_VmaxMif                            
c     par( 38)  :         global_VmaxMMf                            
c     par( 39)  :         global_KiaMi                              
c     par( 40)  :         global_KbMi                               
c     par( 41)  :         global_KicMi                              
c     par( 42)  :         global_KdMi                               
c     par( 43)  :         global_KibMi                              
c     par( 44)  :         global_KidMi                              
c     par( 45)  :         global_KiaMM                              
c     par( 46)  :         global_KbMM                               
c     par( 47)  :         global_KicMM                              
c     par( 48)  :         global_KdMM                               
c     par( 49)  :         global_KibMM                              
c     par( 50)  :         global_KidMM                              
c     par( 51)  :         global_Vmaxsyn                            
c     par( 52)  :         global_Kadp                               
c     par( 53)  :         global_Kpi                                
c     par( 54)  :         global_PSmomPi                            
c     par( 55)  :         global_PSmomCr                            
c     par( 56)  :         global_PSmomPCr                           
c     par( 57)  :         global_pulsatility                        
c     ---
c
c-----------------------------------------------------------------------
c
      double precision    com(3)
      double precision    spe(10)
      double precision    par(57)
      double precision    rul(20)
      double precision    rea(9)
      logical             eve(3)
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
      com(  1) = 6.250000D-02
      com(  2) = 7.500000D-01
      com(  3) = 1.000000D+00
c 
c
      spe(  1) = 3.900000D+01
      spe(  2) = 5.626000D+03
      spe(  3) = 9.789000D+03
      spe(  4) = 5.711000D+03
      spe(  5) = 5.710000D+03
      spe(  6) = 6.400000D+01
      spe(  7) = 5.601000D+03
      spe(  8) = 9.789000D+03
      spe(  9) = 9.100000D+02
      spe( 10) = 9.120000D+02
c 
c
      par(  1) = 1.000000D+00
      par(  2) = 1.000000D+00
      par(  3) = 0.000000D+00
      par(  4) = 0.000000D+00
      par(  5) = 0.000000D+00
      par(  6) = 0.000000D+00
      par(  7) = 1.000000D+00
      par(  8) = 0.000000D+00
      par(  9) = 0.000000D+00
      par( 10) = 0.000000D+00
      par( 11) = 0.000000D+00
      par( 12) = 0.000000D+00
      par( 13) = 1.000000D-03
      par( 14) = 0.000000D+00
      par( 15) = 1.350000D+02
      par( 16) = 1.350000D+02
      par( 17) = 2.200000D+02
      par( 18) = 6.666000D-01
      par( 19) = 1.666667D-01
      par( 20) = 1.666667D-01
      par( 21) = 2.918416D+03
      par( 22) = 3.764847D+03
      par( 23) = 4.000000D+01
      par( 24) = 6.276000D+02
      par( 25) = 4.865000D+02
      par( 26) = 0.000000D+00
      par( 27) = 4.865000D+02
      par( 28) = 1.000000D+00
      par( 29) = 2.860000D-02
      par( 30) = 0.000000D+00
      par( 31) = 1.329477D+01
      par( 32) = 1.520000D+02
      par( 33) = 4.630354D+04
      par( 34) = 3.520341D+03
      par( 35) = 8.820756D+02
      par( 36) = 1.144178D+04
      par( 37) = 8.820756D+02
      par( 38) = 1.144178D+04
      par( 39) = 7.500000D+02
      par( 40) = 5.200000D+03
      par( 41) = 2.048000D+02
      par( 42) = 5.000000D+02
      par( 43) = 2.880000D+04
      par( 44) = 1.600000D+03
      par( 45) = 9.000000D+02
      par( 46) = 1.550000D+04
      par( 47) = 2.224000D+02
      par( 48) = 1.670000D+03
      par( 49) = 3.490000D+04
      par( 50) = 4.730000D+03
      par( 51) = 1.503740D+03
      par( 52) = 2.500000D+01
      par( 53) = 8.000000D+02
      par( 54) = 1.940000D+02
      par( 55) = 1.550000D+02
      par( 56) = 1.550000D+02
      par( 57) = 1.000000D+00
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
      ms = min0(10,ls)
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
         ls = 10
c
      end if
c
c-----------------------------------------------------------------------
c
      mp = min0(57,lp)
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
         lp = 57
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
      idc(  1) = 'IMS'//char(0)
      idc(  2) = 'CYT'//char(0)
      idc(  3) = 'cell'//char(0)
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
      character           ids(10)*(*)
c
      integer             j, nspe
c
c-----------------------------------------------------------------------
c
      nspe = 10
c     ---
      do j = 1, nspe
c
        ids(j) = char(0)
c
      end do
c     ---
      ids(  1) = 'ADPi'//char(0)
      ids(  2) = 'ATPi'//char(0)
      ids(  3) = 'Cri'//char(0)
      ids(  4) = 'PCri'//char(0)
      ids(  5) = 'PCr'//char(0)
      ids(  6) = 'ADP'//char(0)
      ids(  7) = 'ATP'//char(0)
      ids(  8) = 'Cr'//char(0)
      ids(  9) = 'P_ii'//char(0)
      ids( 10) = 'P_i'//char(0)
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
      character           idp(57)*(*)
c
      integer             npar
c
c-----------------------------------------------------------------------
c
      npar = 57
c     ---
      idp(  1) = 'global_j_diff_pcr'//char(0)
      idp(  2) = 'global_j_diff_atp'//char(0)
      idp(  3) = 'global_densyn'//char(0)
      idp(  4) = 'global_tmito'//char(0)
      idp(  5) = 'global_vatpnorm'//char(0)
      idp(  6) = 'global_jsyn'//char(0)
      idp(  7) = 'global_r_diff_pcr'//char(0)
      idp(  8) = 'global_j_ck_mi'//char(0)
      idp(  9) = 'global_j_ck_mm'//char(0)
      idp( 10) = 'global_j_diff_adp'//char(0)
      idp( 11) = 'global_j_diff_cr'//char(0)
      idp( 12) = 'global_j_diff_pi'//char(0)
      idp( 13) = 'global_stepsize'//char(0)
      idp( 14) = 'global_phase'//char(0)
      idp( 15) = 'global_heartrate_bpm'//char(0)
      idp( 16) = 'global_heartrate_basis'//char(0)
      idp( 17) = 'global_heartrate_test'//char(0)
      idp( 18) = 'global_fracDia'//char(0)
      idp( 19) = 'global_fracSysUp'//char(0)
      idp( 20) = 'global_fracSysDown'//char(0)
      idp( 21) = 'global_VhydAmp_basis'//char(0)
      idp( 22) = 'global_VhydAmp_test'//char(0)
      idp( 23) = 'global_time_Jhyd_step'//char(0)
      idp( 24) = 'global_Jhyd_test'//char(0)
      idp( 25) = 'global_Jhyd_basis'//char(0)
      idp( 26) = 'global_last_time_fired'//char(0)
      idp( 27) = 'global_Jhyd'//char(0)
      idp( 28) = 'global_ck_factor_iaa'//char(0)
      idp( 29) = 'global_ck_factor_ia'//char(0)
      idp( 30) = 'global_tmito_factor'//char(0)
      idp( 31) = 'global_PSmomATP'//char(0)
      idp( 32) = 'global_K_CK_eq'//char(0)
      idp( 33) = 'global_VmaxMMb'//char(0)
      idp( 34) = 'global_VmaxMib'//char(0)
      idp( 35) = 'global_VmaxMif_full_activity'//char(0)
      idp( 36) = 'global_VmaxMMf_full_activity'//char(0)
      idp( 37) = 'global_VmaxMif'//char(0)
      idp( 38) = 'global_VmaxMMf'//char(0)
      idp( 39) = 'global_KiaMi'//char(0)
      idp( 40) = 'global_KbMi'//char(0)
      idp( 41) = 'global_KicMi'//char(0)
      idp( 42) = 'global_KdMi'//char(0)
      idp( 43) = 'global_KibMi'//char(0)
      idp( 44) = 'global_KidMi'//char(0)
      idp( 45) = 'global_KiaMM'//char(0)
      idp( 46) = 'global_KbMM'//char(0)
      idp( 47) = 'global_KicMM'//char(0)
      idp( 48) = 'global_KdMM'//char(0)
      idp( 49) = 'global_KibMM'//char(0)
      idp( 50) = 'global_KidMM'//char(0)
      idp( 51) = 'global_Vmaxsyn'//char(0)
      idp( 52) = 'global_Kadp'//char(0)
      idp( 53) = 'global_Kpi'//char(0)
      idp( 54) = 'global_PSmomPi'//char(0)
      idp( 55) = 'global_PSmomCr'//char(0)
      idp( 56) = 'global_PSmomPCr'//char(0)
      idp( 57) = 'global_pulsatility'//char(0)
c     ---
c
c-----------------------------------------------------------------------
c
      return
      end subroutine get_parameter_ids
c
c=======================================================================
c
c 
c
c=======================================================================
c
      double precision function piecewise2 ( val1, cond1, dflt )
c 
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    val1
      double precision    dflt
c 
      logical             cond1
c 
c-----------------------------------------------------------------------
c
      if ( cond1 ) then
c
         piecewise2 = val1
c 
c 
      else
c
         piecewise2 = dflt
c
      end if
c
c-----------------------------------------------------------------------
c
      return
      end function piecewise2
c
c=======================================================================
c
      double precision function piecewise6 ( val1, cond1, val2,
     &          cond2, val3, cond3, dflt )
c 
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    val1
      double precision    val2
      double precision    val3
      double precision    dflt
c 
      logical             cond1
      logical             cond2
      logical             cond3
c 
c-----------------------------------------------------------------------
c
      if ( cond1 ) then
c
         piecewise6 = val1
c 
      else if ( cond2 ) then
c
         piecewise6 = val2
c 
      else if ( cond3 ) then
c
         piecewise6 = val3
c 
      else
c
         piecewise6 = dflt
c
      end if
c
c-----------------------------------------------------------------------
c
      return
      end function piecewise6
c
c 
c

