#ifndef __SBML_MATH_CPP_HPP
#define __SBML_MATH_CPP_HPP

#define NDEBUG
#include <cmath>
#include <cassert>

namespace ODEydot
{
  template<typename T>   T         abs( T x );
  template<typename T>   T        acos( T x );
  template<typename T>   T      arccos( T x );
  template<typename T>   T       acosh( T x );
  template<typename T>   T     arccosh( T x );
  template<typename T>   T        acot( T x );
  template<typename T>   T      arccot( T x );
  template<typename T>   T       acoth( T x );
  template<typename T>   T     arccoth( T x );
  template<typename T>   T        acsc( T x );
  template<typename T>   T      arccsc( T x );
  template<typename T>   T       acsch( T x );
  template<typename T>   T     arccsch( T x );
  template<typename T>   T        asec( T x );
  template<typename T>   T      arcsec( T x );
  template<typename T>   T       asech( T x );
  template<typename T>   T     arcsech( T x );
  template<typename T>   T        asin( T x );
  template<typename T>   T      arcsin( T x );
  template<typename T>   T       asinh( T x );
  template<typename T>   T     arcsinh( T x );
  template<typename T>   T        atan( T x );
  template<typename T>   T      arctan( T x );
  template<typename T>   T       atanh( T x );
  template<typename T>   T     arctanh( T x );

  template<typename T>   T        ceil( T x );
  template<typename T>   T     ceiling( T x );
  template<typename T>   T         cos( T x );
  template<typename T>   T        cosh( T x );
  template<typename T>   T         cot( T x );
  template<typename T>   T        coth( T x );
  template<typename T>   T         csc( T x );
  template<typename T>   T        csch( T x );

  template<typename T>   T         exp( T x );

  template<typename T>   T   factorial( int x );
  template<typename T>   T       floor( T x );

  template<typename T>   T          ln( T x );
  template<typename T>   T         log( T x );
  template<typename T>   T         log( T x, T y );
  template<typename T>   T       log10( T x );

  template<typename T>   T         pow( T x, T y ); 
  template<typename T>   T       power( T x, T y ); 

  template<typename T>   T        root( T x, T y ); 

  template<typename T>   T         sec( T x );
  template<typename T>   T        sech( T x );
  template<typename T>   T         sqr( T x );
  template<typename T>   T        sqrt( T x );
  template<typename T>   T         sin( T x );
  template<typename T>   T        sinh( T x );

  template<typename T>   T         tan( T x );
  template<typename T>   T        tanh( T x );
};

//=======================================================================
template<typename T>
T ODEydot::abs( T x )
{
  return ::fabs(x);
}
//=======================================================================
template<typename T>
T ODEydot::exp( T x )
{
  return ::exp(x);
}
//=======================================================================
template<typename T>
T ODEydot::log( T x )
{
  return ::log(x);
}
//=======================================================================
template<typename T>
T ODEydot::pow( T x, T y )
{
  assert( x >= 0.0 );
  T ret = ::pow(::fabs(x),y);
  // T ret = ::exp( y*::log(::fabs(x)) );
  return ret;
}
//=======================================================================

#endif // __SBML_MATH_CPP_HPP

/*
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
*/
