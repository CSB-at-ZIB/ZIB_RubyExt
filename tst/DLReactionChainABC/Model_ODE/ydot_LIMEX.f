c-----
c SBML Model : ChainABC                                               
c              ~~~~~~~~
c       Date : Sun Apr  5 16:58:44 2015
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
c     spe(  1)  :         A                                         
c     spe(  2)  :         B                                         
c     spe(  3)  :         C                                         
c     ---
c
c     ---
c     par(  1)  :         global_k1                                 
c     par(  2)  :         global_k_1                                
c     par(  3)  :         global_k2                                 
c     par(  4)  :         global_k_2                                
c     ---
c
      double precision    com(1)
      double precision    spe(3)
      double precision    par(4)
      double precision    rul(0)
      double precision    rea(4)
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
c 
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      rea(1) = par(1) * spe(1)
      rea(2) = par(2) * spe(2)
      rea(3) = par(3) * spe(2)
      rea(4) = par(4) * spe(3)
c 
c-----------------------------------------------------------------------
c
      call check_events(t, 0, eve)
c
c-----------------------------------------------------------------------
c
      dy(1) = ( - rea(1) + rea(2) ) / com(1)
      dy(2) = ( + rea(1) - rea(2) - rea(3) + rea(4) ) / com(1)
      dy(3) = ( + rea(3) - rea(4) ) / com(1)
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
      double precision    spe(3)
      double precision    par(4)
      double precision    rul(0)
      double precision    rea(4)
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
c     ---
c
c     ---
c     spe(  1)  :         A                                         
c     spe(  2)  :         B                                         
c     spe(  3)  :         C                                         
c     ---
c
c     ---
c     par(  1)  :         global_k1                                 
c     par(  2)  :         global_k_1                                
c     par(  3)  :         global_k2                                 
c     par(  4)  :         global_k_2                                
c     ---
c
c-----------------------------------------------------------------------
c
      double precision    com(1)
      double precision    spe(3)
      double precision    par(4)
      double precision    rul(0)
      double precision    rea(4)
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
      spe(  1) = 1.000000D+00
      spe(  2) = 0.000000D+00
      spe(  3) = 0.000000D+00
c 
c
      par(  1) = 2.000000D+00
      par(  2) = 3.000000D-03
      par(  3) = 1.000000D+00
      par(  4) = 2.000000D-03
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
      ms = min0(3,ls)
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
         ls = 3
c
      end if
c
c-----------------------------------------------------------------------
c
      mp = min0(4,lp)
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
         lp = 4
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
      character           ids(3)*(*)
c
      integer             j, nspe
c
c-----------------------------------------------------------------------
c
      nspe = 3
c     ---
      do j = 1, nspe
c
        ids(j) = char(0)
c
      end do
c     ---
      ids(  1) = 'A'//char(0)
      ids(  2) = 'B'//char(0)
      ids(  3) = 'C'//char(0)
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
      character           idp(4)*(*)
c
      integer             j, npar
c
c-----------------------------------------------------------------------
c
      npar = 4
c     ---
      do j = 1, npar
c
        idp(j) = char(0)
c
      end do
c     ---
      idp(  1) = 'global_k1'//char(0)
      idp(  2) = 'global_k_1'//char(0)
      idp(  3) = 'global_k2'//char(0)
      idp(  4) = 'global_k_2'//char(0)
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
      idm(  1) = 'ChainABC'//char(0)
      idm(  2) = 'Sun Apr  5 16:58:44 2015'//char(0)
      idm(  3) = '001428245924'//char(0)
      idm(  4) = 'ChainABC.xml'//char(0)
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

