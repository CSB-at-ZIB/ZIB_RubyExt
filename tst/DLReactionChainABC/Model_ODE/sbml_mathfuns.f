c
c=======================================================================
c
      double precision function abs ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      abs = dabs(x)
c
c-----------------------------------------------------------------------
c
      return
      end function abs
c
c=======================================================================
c
      double precision function acos ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      acos = dacos(x)
c
c-----------------------------------------------------------------------
c
      return
      end function acos
c
c=======================================================================
c
      double precision function arccos ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      arccos = dacos(x)
c
c-----------------------------------------------------------------------
c
      return
      end function arccos
c
c=======================================================================
c
      double precision function acosh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      acosh = dacosh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function acosh
c
c=======================================================================
c
      double precision function arccosh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      arccosh = dacosh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function arccosh
c
c=======================================================================
c
      double precision function acot ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one, pihalf
c
      parameter           ( one = 1.0d0 , pihalf = 2.0*datan(one) )
c
c-----------------------------------------------------------------------
c
      acot = pihalf - datan(x)
c
c-----------------------------------------------------------------------
c
      return
      end function acot
c
c=======================================================================
c
      double precision function arccot ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one, pihalf
c
      parameter           ( one = 1.0d0 , pihalf = 2.0*datan(one) )
c
c-----------------------------------------------------------------------
c
      arccot = pihalf - datan(x)
c
c-----------------------------------------------------------------------
c
      return
      end function arccot
c
c=======================================================================
c
      double precision function acoth ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      acoth = datanh(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function acoth
c
c=======================================================================
c
      double precision function arccoth ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      arccoth = datanh(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function arccoth
c
c=======================================================================
c
      double precision function acsc ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      acsc = dasin(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function acsc
c
c=======================================================================
c
      double precision function arccsc ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      arccsc = dasin(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function arccsc
c
c=======================================================================
c
      double precision function acsch ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      acsch = dasinh(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function acsch
c
c=======================================================================
c
      double precision function arccsch ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      arccsch = dasinh(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function arccsch
c
c=======================================================================
c
      double precision function asec ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      asec = dacos(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function asec
c
c=======================================================================
c
      double precision function arcsec ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      arcsec = dacos(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function arcsec
c
c=======================================================================
c
      double precision function asech ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      asech = dacosh(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function asech
c
c=======================================================================
c
      double precision function arcsech ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      arcsech = dacosh(one / x)
c
c-----------------------------------------------------------------------
c
      return
      end function arcsech
c
c=======================================================================
c
      double precision function asin ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      asin = dasin(x)
c
c-----------------------------------------------------------------------
c
      return
      end function asin
c
c=======================================================================
c
      double precision function arcsin ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      arcsin = dasin(x)
c
c-----------------------------------------------------------------------
c
      return
      end function arcsin
c
c=======================================================================
c
      double precision function asinh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      asinh = dasinh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function asinh
c
c=======================================================================
c
      double precision function arcsinh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      arcsinh = dasinh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function arcsinh
c
c=======================================================================
c
      double precision function atan ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      atan = datan(x)
c
c-----------------------------------------------------------------------
c
      return
      end function atan
c
c=======================================================================
c
      double precision function arctan ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      arctan = datan(x)
c
c-----------------------------------------------------------------------
c
      return
      end function arctan
c
c=======================================================================
c
      double precision function atanh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      atanh = datanh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function atanh
c
c=======================================================================
c
      double precision function arctanh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      arctanh = datanh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function arctanh
c
c=======================================================================
c
      double precision function ceil ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      ceil = dnint(x)
c
c-----------------------------------------------------------------------
c
      return
      end function ceil
c
c=======================================================================
c
      double precision function ceiling ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      ceiling = dnint(x)
c
c-----------------------------------------------------------------------
c
      return
      end function ceiling
c
c=======================================================================
c
      double precision function cos ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      cos = dcos(x)
c
c-----------------------------------------------------------------------
c
      return
      end function cos
c
c=======================================================================
c
      double precision function cosh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      cosh = dcosh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function cosh
c
c=======================================================================
c
      double precision function cot ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      cot = dcos(x) / dsin(x)
c
c-----------------------------------------------------------------------
c
      return
      end function cot
c
c=======================================================================
c
      double precision function coth ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      coth = dcosh(x) / dsinh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function coth
c
c=======================================================================
c
      double precision function csc ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      csc = one / dsin(x)
c
c-----------------------------------------------------------------------
c
      return
      end function csc
c
c=======================================================================
c
      double precision function csch ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 ) 
c
c-----------------------------------------------------------------------
c
      csch = one / dsinh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function csch
c
c=======================================================================
c
      double precision function factorial ( n )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      integer             n, j
c
      double precision    one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      factorial = one
c
      do j = 2, n
c
         factorial = factorial*j
c
      end do
c
c-----------------------------------------------------------------------
c
      return
      end function factorial
c
c=======================================================================
c
      double precision function exp ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      exp = dexp(x)
c
c-----------------------------------------------------------------------
c
      return
      end function exp
c
c=======================================================================
c
      double precision function floor ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      floor = dint(x)
c
c-----------------------------------------------------------------------
c
      return
      end function floor
c
c=======================================================================
c
      double precision function ln ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      ln = dlog(x)
c
c-----------------------------------------------------------------------
c
      return
      end function ln
c
c=======================================================================
c
      double precision function log ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, y
c
c-----------------------------------------------------------------------
c
      log = dlog(y) / dlog(x)
c
c-----------------------------------------------------------------------
c
      return
      end function log
c
c=======================================================================
c
      double precision function log10 ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, ten, logten
c
      parameter           ( ten = 1.0d1 , logten = dlog(ten) )
c
c-----------------------------------------------------------------------
c
      log10 = dlog(x) / logten
c
c-----------------------------------------------------------------------
c
      return
      end function log10
c
c=======================================================================
c
      double precision function pow ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, y
c
c-----------------------------------------------------------------------
c
      pow = x**y
c
c-----------------------------------------------------------------------
c
      return
      end function pow
c
c=======================================================================
c
      double precision function power ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, y
c
c-----------------------------------------------------------------------
c
      power = x**y
c
c-----------------------------------------------------------------------
c
      return
      end function power
c
c=======================================================================
c
      double precision function root ( b, x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    b, x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      root = x**(one/b)
c
c-----------------------------------------------------------------------
c
      return
      end function root
c
c=======================================================================
c
      double precision function sec ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      sec = one / dcos(x)
c
c-----------------------------------------------------------------------
c
      return
      end function sec
c
c=======================================================================
c
      double precision function sech ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x, one
c
      parameter           ( one = 1.0d0 )
c
c-----------------------------------------------------------------------
c
      sech = one / dcosh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function sech
c
c=======================================================================
c
      double precision function sqr ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      sqr = x*x
c
c-----------------------------------------------------------------------
c
      return
      end function sqr
c
c=======================================================================
c
      double precision function sqrt ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      sqrt = dsqrt(x)
c
c-----------------------------------------------------------------------
c
      return
      end function sqrt
c
c=======================================================================
c
      double precision function sin ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      sin = dsin(x)
c
c-----------------------------------------------------------------------
c
      return
      end function sin
c
c=======================================================================
c
      double precision function sinh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      sinh = dsinh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function sinh
c
c=======================================================================
c
      double precision function tan ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      tan = dtan(x)
c
c-----------------------------------------------------------------------
c
      return
      end function tan
c
c=======================================================================
c
      double precision function tanh ( x )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
c
c-----------------------------------------------------------------------
c
      tanh = dtanh(x)
c
c-----------------------------------------------------------------------
c
      return
      end function tanh
c
c=======================================================================
c
c=======================================================================
c
      logical function and ( a, b )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      logical             a
      logical             b
c
c-----------------------------------------------------------------------
c
      and = a.AND.b
c
c-----------------------------------------------------------------------
c
      return
      end function and
c
c=======================================================================
c
      logical function or ( a, b )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      logical             a
      logical             b
c
c-----------------------------------------------------------------------
c
      or = a.OR.b
c
c-----------------------------------------------------------------------
c
      return
      end function or
c
c=======================================================================
c
      logical function not ( a )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      logical             a
c
c-----------------------------------------------------------------------
c
      not = .NOT.a
c
c-----------------------------------------------------------------------
c
      return
      end function not
c
c=======================================================================
c 
      logical function eq ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
      double precision    y
c
c-----------------------------------------------------------------------
c
      eq = x.EQ.y
c
c-----------------------------------------------------------------------
c
      return
      end function eq
c
c=======================================================================
c 
      logical function neq ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
      double precision    y
c
c-----------------------------------------------------------------------
c
      neq = x.NE.y
c
c-----------------------------------------------------------------------
c
      return
      end function neq
c
c=======================================================================
c 
      logical function gt ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
      double precision    y
c
c-----------------------------------------------------------------------
c
      gt = x.GT.y
c
c-----------------------------------------------------------------------
c
      return
      end function gt
c
c=======================================================================
c 
      logical function geq ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
      double precision    y
c
c-----------------------------------------------------------------------
c
      geq = x.GE.y
c
c-----------------------------------------------------------------------
c
      return
      end function geq
c
c=======================================================================
c 
      logical function lt ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
      double precision    y
c
c-----------------------------------------------------------------------
c
      lt = x.LT.y
c
c-----------------------------------------------------------------------
c
      return
      end function lt
c
c=======================================================================
c 
      logical function leq ( x, y )
c
      implicit none
c
c-----------------------------------------------------------------------
c
      double precision    x
      double precision    y
c
c-----------------------------------------------------------------------
c
      leq = x.LE.y
c
c-----------------------------------------------------------------------
c
      return
      end function leq
c
c=======================================================================
c 
