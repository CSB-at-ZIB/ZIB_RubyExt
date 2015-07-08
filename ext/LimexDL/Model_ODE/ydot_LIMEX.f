c-----
c SBML Model : Gardner1998_CellCycle_Goldbeter                        
c              ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
c       Date : Wed Jan 21 18:33:42 2015
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
c     com(  2)  :         Cell                                      
c     ---
c
c     ---
c     spe(  1)  :         C                                         
c     spe(  2)  :         X                                         
c     spe(  3)  :         M                                         
c     spe(  4)  :         Y                                         
c     spe(  5)  :         Z                                         
c     ---
c
c     ---
c     par(  1)  :         global_V1                                 
c     par(  2)  :         global_K6                                 
c     par(  3)  :         global_V1p                                
c     par(  4)  :         global_V3                                 
c     par(  5)  :         global_V3p                                
c     par(  6)  :         reaction1_vi                              
c     par(  7)  :         reaction2_k1                              
c     par(  8)  :         reaction2_K5                              
c     par(  9)  :         reaction3_kd                              
c     par( 10)  :         reaction4_K1                              
c     par( 11)  :         reaction5_V2                              
c     par( 12)  :         reaction5_K2                              
c     par( 13)  :         reaction6_K3                              
c     par( 14)  :         reaction7_K4                              
c     par( 15)  :         reaction7_V4                              
c     par( 16)  :         reaction8_a1                              
c     par( 17)  :         reaction9_a2                              
c     par( 18)  :         reaction10_alpha                          
c     par( 19)  :         reaction10_d1                             
c     par( 20)  :         reaction11_kd                             
c     par( 21)  :         reaction11_alpha                          
c     par( 22)  :         reaction12_vs                             
c     par( 23)  :         reaction13_d1                             
c     ---
c
      double precision    com(2)
      double precision    spe(5)
      double precision    par(23)
      double precision    rul(2)
      double precision    rea(13)
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
      rul(1) = spe(1) * par(3) * pow(spe(1) + par(2),
     &          (-1.000000D+00))
      rul(2) = spe(3) * par(5)
c 
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      rea(1) = par(6)
      rea(2) = spe(1) * par(7) * spe(2) * pow(spe(1) + par(8),
     &          (-1.000000D+00))
      rea(3) = spe(1) * par(9)
      rea(4) = (1.000000D+00 + (-1.000000D+00) * spe(3)) * rul(1) *
     &          pow(par(10) + (-1.000000D+00) * spe(3) +
     &          1.000000D+00, (-1.000000D+00))
      rea(5) = spe(3) * par(11) * pow(par(12) + spe(3),
     &          (-1.000000D+00))
      rea(6) = rul(2) * (1.000000D+00 + (-1.000000D+00) * spe(2)) *
     &          pow(par(13) + (-1.000000D+00) * spe(2) +
     &          1.000000D+00, (-1.000000D+00))
      rea(7) = par(15) * spe(2) * pow(par(14) + spe(2),
     &          (-1.000000D+00))
      rea(8) = par(16) * spe(1) * spe(4)
      rea(9) = par(17) * spe(5)
      rea(10) = par(18) * par(19) * spe(5)
      rea(11) = par(21) * par(20) * spe(5)
      rea(12) = par(22)
      rea(13) = par(23) * spe(4)
c 
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      dy(1) = ( + rea(1) - rea(2) - rea(3) - rea(8) + rea(9) +
     &          rea(10) ) / com(2)
      dy(2) = ( + rea(6) - rea(7) ) / com(2)
      dy(3) = ( + rea(4) - rea(5) ) / com(2)
      dy(4) = ( - rea(8) + rea(9) + rea(11) + rea(12) - rea(13) ) /
     &          com(2)
      dy(5) = ( + rea(8) - rea(9) - rea(10) - rea(11) ) / com(2)
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
      double precision    com(2)
      double precision    spe(5)
      double precision    par(23)
      double precision    rul(2)
      double precision    rea(13)
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
c     com(  2)  :         Cell                                      
c     ---
c
c     ---
c     spe(  1)  :         C                                         
c     spe(  2)  :         X                                         
c     spe(  3)  :         M                                         
c     spe(  4)  :         Y                                         
c     spe(  5)  :         Z                                         
c     ---
c
c     ---
c     par(  1)  :         global_V1                                 
c     par(  2)  :         global_K6                                 
c     par(  3)  :         global_V1p                                
c     par(  4)  :         global_V3                                 
c     par(  5)  :         global_V3p                                
c     par(  6)  :         reaction1_vi                              
c     par(  7)  :         reaction2_k1                              
c     par(  8)  :         reaction2_K5                              
c     par(  9)  :         reaction3_kd                              
c     par( 10)  :         reaction4_K1                              
c     par( 11)  :         reaction5_V2                              
c     par( 12)  :         reaction5_K2                              
c     par( 13)  :         reaction6_K3                              
c     par( 14)  :         reaction7_K4                              
c     par( 15)  :         reaction7_V4                              
c     par( 16)  :         reaction8_a1                              
c     par( 17)  :         reaction9_a2                              
c     par( 18)  :         reaction10_alpha                          
c     par( 19)  :         reaction10_d1                             
c     par( 20)  :         reaction11_kd                             
c     par( 21)  :         reaction11_alpha                          
c     par( 22)  :         reaction12_vs                             
c     par( 23)  :         reaction13_d1                             
c     ---
c
c-----------------------------------------------------------------------
c
      double precision    com(2)
      double precision    spe(5)
      double precision    par(23)
      double precision    rul(2)
      double precision    rea(13)
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
      com(  2) = 1.000000D+00
c 
c
      spe(  1) = 0.000000D+00
      spe(  2) = 0.000000D+00
      spe(  3) = 0.000000D+00
      spe(  4) = 1.000000D+00
      spe(  5) = 1.000000D+00
c 
c
      par(  1) = 0.000000D+00
      par(  2) = 3.000000D-01
      par(  3) = 7.500000D-01
      par(  4) = 0.000000D+00
      par(  5) = 3.000000D-01
      par(  6) = 1.000000D-01
      par(  7) = 5.000000D-01
      par(  8) = 5.000000D-02
      par(  9) = 2.000000D-02
      par( 10) = 1.000000D-01
      par( 11) = 2.500000D-01
      par( 12) = 5.000000D-01
      par( 13) = 5.000000D-01
      par( 14) = 5.000000D-01
      par( 15) = 1.000000D-01
      par( 16) = 5.000000D-02
      par( 17) = 5.000000D-02
      par( 18) = 1.000000D-01
      par( 19) = 5.000000D-02
      par( 20) = 2.000000D-02
      par( 21) = 1.000000D-01
      par( 22) = 2.000000D-01
      par( 23) = 5.000000D-02
c 
c
c 
c
c-----------------------------------------------------------------------
c
      mc = min0(2,lc)
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
         lc = 2
c
      end if
c
c-----------------------------------------------------------------------
c
      ms = min0(5,ls)
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
         ls = 5
c
      end if
c
c-----------------------------------------------------------------------
c
      mp = min0(23,lp)
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
         lp = 23
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
      character           idc(2)*(*)
c
      integer             ncom
c
c-----------------------------------------------------------------------
c
      ncom = 2
c     ---
      idc(  1) = 'default'//char(0)
      idc(  2) = 'Cell'//char(0)
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
      character           ids(5)*(*)
c
      integer             j, nspe
c
c-----------------------------------------------------------------------
c
      nspe = 5
c     ---
      do j = 1, nspe
c
        ids(j) = char(0)
c
      end do
c     ---
      ids(  1) = 'C'//char(0)
      ids(  2) = 'X'//char(0)
      ids(  3) = 'M'//char(0)
      ids(  4) = 'Y'//char(0)
      ids(  5) = 'Z'//char(0)
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
      character           idp(23)*(*)
c
      integer             j, npar
c
c-----------------------------------------------------------------------
c
      npar = 23
c     ---
      do j = 1, npar
c
        idp(j) = char(0)
c
      end do
c     ---
      idp(  1) = 'global_V1'//char(0)
      idp(  2) = 'global_K6'//char(0)
      idp(  3) = 'global_V1p'//char(0)
      idp(  4) = 'global_V3'//char(0)
      idp(  5) = 'global_V3p'//char(0)
      idp(  6) = 'reaction1_vi'//char(0)
      idp(  7) = 'reaction2_k1'//char(0)
      idp(  8) = 'reaction2_K5'//char(0)
      idp(  9) = 'reaction3_kd'//char(0)
      idp( 10) = 'reaction4_K1'//char(0)
      idp( 11) = 'reaction5_V2'//char(0)
      idp( 12) = 'reaction5_K2'//char(0)
      idp( 13) = 'reaction6_K3'//char(0)
      idp( 14) = 'reaction7_K4'//char(0)
      idp( 15) = 'reaction7_V4'//char(0)
      idp( 16) = 'reaction8_a1'//char(0)
      idp( 17) = 'reaction9_a2'//char(0)
      idp( 18) = 'reaction10_alpha'//char(0)
      idp( 19) = 'reaction10_d1'//char(0)
      idp( 20) = 'reaction11_kd'//char(0)
      idp( 21) = 'reaction11_alpha'//char(0)
      idp( 22) = 'reaction12_vs'//char(0)
      idp( 23) = 'reaction13_d1'//char(0)
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
      character           idm(5)*(*)
c
      integer             j, nid
c
c-----------------------------------------------------------------------
c
      nid = 5
c     ---
      do j = 1, nid
c
        idm(j) = char(0)
c
      end do
c     ---
      idm(  1) = 'Gardner1998_CellCycle_Goldbeter'//char(0)
      idm(  2) = 'Wed Jan 21 18:33:42 2015'//char(0)
      idm(  3) = '001421861622'//char(0)
      idm(  4) = 'BIOMD0000000008.xml'//char(0)
      idm(  5) = 'no var eq'//char(0)
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

