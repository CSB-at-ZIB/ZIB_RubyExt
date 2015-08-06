#ifndef __SBML_MATH_CPP_HPP
#define __SBML_MATH_CPP_HPP

#define NDEBUG
#include <cassert>
#include <cmath>
//=======================================================================
namespace ODEydot
{
  template<typename T>   T        abs( T x );
  template<typename T>   T       acos( T x );
  template<typename T>   T     arccos( T x );
  template<typename T>   T      acosh( T x );
  template<typename T>   T    arccosh( T x );
  template<typename T>   T       acot( T x );
  template<typename T>   T     arccot( T x );
  template<typename T>   T      acoth( T x );
  template<typename T>   T    arccoth( T x );
  template<typename T>   T       acsc( T x );
  template<typename T>   T     arccsc( T x );
  template<typename T>   T      acsch( T x );
  template<typename T>   T    arccsch( T x );
  template<typename T>   T       asec( T x );
  template<typename T>   T     arcsec( T x );
  template<typename T>   T      asech( T x );
  template<typename T>   T    arcsech( T x );
  template<typename T>   T       asin( T x );
  template<typename T>   T     arcsin( T x );
  template<typename T>   T      asinh( T x );
  template<typename T>   T    arcsinh( T x );
  template<typename T>   T       atan( T x );
  template<typename T>   T     arctan( T x );
  template<typename T>   T      atanh( T x );
  template<typename T>   T    arctanh( T x );

  template<typename T>   T       ceil( T x );
  template<typename T>   T    ceiling( T x );
  template<typename T>   T        cos( T x );
  template<typename T>   T       cosh( T x );
  template<typename T>   T        cot( T x );
  template<typename T>   T       coth( T x );
  template<typename T>   T        csc( T x );
  template<typename T>   T       csch( T x );

  template<typename T>   T        exp( T x );

  template<typename T>   T  factorial( T x );
  template<typename T>   T      floor( T x );

  template<typename T>   T         ln( T x );
  template<typename T>   T        log( T x );
  template<typename T>   T        log( T b, T x );
  template<typename T>   T      log10( T x );

  template<typename T>   T        pow( T x, T y ); 
  template<typename T>   T      power( T x, T y ); 

  template<typename T>   T       root( T b, T x ); 

  template<typename T>   T        sec( T x );
  template<typename T>   T       sech( T x );
  template<typename T>   T        sqr( T x );
  template<typename T>   T       sqrt( T x );
  template<typename T>   T        sin( T x );
  template<typename T>   T       sinh( T x );

  template<typename T>   T        tan( T x );
  template<typename T>   T       tanh( T x );

  template<typename T>   int sbml_and( T a, T b );
  template<typename T>   int sbml_or ( T a, T b );
  template<typename T>   int sbml_not( T a );

  template<typename T1, typename T2>   int    eq( T1 x, T2 y );
  template<typename T1, typename T2>   int   neq( T1 x, T2 y );
  template<typename T1, typename T2>   int    gt( T1 x, T2 y );
  template<typename T1, typename T2>   int   geq( T1 x, T2 y );
  template<typename T1, typename T2>   int    lt( T1 x, T2 y );
  template<typename T1, typename T2>   int   leq( T1 x, T2 y );
};
//=======================================================================
template<typename T>
T ODEydot::abs( T x )
{
  return ::fabs(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::acos( T x )
{
  return ::acos(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arccos( T x )
{
  return ::acos(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::acosh( T x )
{
  return ::acosh(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arccosh( T x )
{
  return ::acosh(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::acot( T x )
{
  return (M_PI_2 - ::atan(x));  // M_PI_2 = M_PI/2
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arccot( T x )
{
  return (M_PI_2 - ::atan(x));  // M_PI_2 = M_PI/2
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::acoth( T x )
{
  assert( x != 0.0 );
  return ::atanh( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arccoth( T x )
{
  assert( x != 0.0 );
  return ::atanh( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::acsc( T x )
{
  assert( x != 0.0 );
  return ::asin( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arccsc( T x )
{
  assert( x != 0.0 );
  return ::asin( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::acsch( T x )
{
  assert( x != 0.0 );
  return ::asinh( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arccsch( T x )
{
  assert( x != 0.0 );
  return ::asinh( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::asec( T x )
{
  assert( x != 0.0 );
  return ::acos( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arcsec( T x )
{
  assert( x != 0.0 );
  return ::acos( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::asech( T x )
{
  assert( x != 0.0 );
  return ::acosh( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arcsech( T x )
{
  assert( x != 0.0 );
  return ::acosh( 1.0 / x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::asin( T x )
{
  return ::asin(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arcsin( T x )
{
  return ::asin(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::asinh( T x )
{
  return ::asinh(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arcsinh( T x )
{
  return ::asinh(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::atan( T x )
{
  return ::atan(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arctan( T x )
{
  return ::atan(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::atanh( T x )
{
  return ::atanh(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::arctanh( T x )
{
  return ::atanh(x);
}
//=======================================================================
template<typename T>
T ODEydot::ceil( T x )
{
  return ::ceil(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::ceiling( T x )
{
  return ::ceil(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::cos( T x )
{
  return ::cos(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::cosh( T x )
{
  return ::cosh(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::cot( T x )
{
  assert( ::sin(x) != 0.0 );
  return ( ::cos(x) / ::sin(x) );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::coth( T x )
{
  assert( x != 0.0 );
  return ( ::cosh(x) / ::sinh(x) );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::csc( T x )
{
  assert( ::sin(x) != 0.0 );
  return ( 1.0 / ::sin(x) );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::csch( T x )
{
  assert( x != 0.0 );
  return ( 1.0 / ::sinh(x) );
}
//=======================================================================
template<typename T>
T ODEydot::factorial( T x )
{
  assert( x > -1.0 );
  return ::gamma( x + 1.0 );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::floor( T x )
{
  return ::floor(x);
}
//=======================================================================
template<typename T>
T ODEydot::exp( T x )
{
  return ::exp(x);
}
//=======================================================================
template<typename T>
T ODEydot::ln( T x )
{
  assert( x > 0.0 );
  return ::log(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::log( T x )
{
  assert( x > 0.0 );
  return ::log(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::log( T b, T x )
{
  assert( x > 0.0 );
  assert( (b > 0.0) && (b != 1.0) );
  return ( ::log(x) / ::log(b) );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::log10( T x )
{
  assert( x > 0.0 );
  return ( ::log(x) / M_LN10 );  // M_LN10 = log(10.0)
}
//=======================================================================
template<typename T>
T ODEydot::pow( T x, T y )
{
  assert( x >= 0.0 );
  T ret = (x != 0.0) ? ::pow( x, y ) : 0.0;
  // T ret = ::exp( y * ::log(::fabs(x)) );
  return ret;
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::power( T x, T y )
{
  assert( x >= 0.0 );
  T ret = (x != 0.0) ? ::pow( x, y ) : 0.0;
  return ret;
}
//=======================================================================
template<typename T>
T ODEydot::root( T b, T x )
{
  assert( x >= 0.0 );
  assert( b != 0.0 );
  T ret = (x != 0.0) ? ::pow( x, 1.0/b ) : 0.0;
  return ret;
}
//=======================================================================
template<typename T>
T ODEydot::sec( T x )
{
  assert( ::cos(x) != 0.0 );
  return ( 1.0 / ::cos(x) );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::sech( T x )
{
  return ( 1.0 / ::cosh(x) );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::sqr( T x )
{
  return ( x*x );
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::sqrt( T x )
{
  return ::sqrt(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::sin( T x )
{
  return ::sin(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::sinh( T x )
{
  return ::sinh(x);
}
//=======================================================================
template<typename T>
T ODEydot::tan( T x )
{
  return ::tan(x);
}
//-----------------------------------------------------------------------
template<typename T>
T ODEydot::tanh( T x )
{
  return ::tanh(x);
}
//=======================================================================
//-----------------------------------------------------------------------
//=======================================================================
template<typename T>
int ODEydot::sbml_and( T a, T b )
{
  return (a and b);
}
//-----------------------------------------------------------------------
template<typename T>
int ODEydot::sbml_or( T a, T b )
{
  return (a or b);
}
//-----------------------------------------------------------------------
template<typename T>
int ODEydot::sbml_not( T a )
{
  return (not a);
}
//=======================================================================
//-----------------------------------------------------------------------
//=======================================================================
template<typename T1, typename T2>
int ODEydot::eq( T1 x, T2 y )
{
  return ( x == y );
}
//-----------------------------------------------------------------------
template<typename T1, typename T2>
int ODEydot::neq( T1 x, T2 y )
{
  return ( x != y );
}
//-----------------------------------------------------------------------
template<typename T1, typename T2>
int ODEydot::gt( T1 x, T2 y )
{
  return ( x > y );
}
//-----------------------------------------------------------------------
template<typename T1, typename T2>
int ODEydot::geq( T1 x, T2 y )
{
  return ( x >= y );
}
//-----------------------------------------------------------------------
template<typename T1, typename T2>
int ODEydot::lt( T1 x, T2 y )
{
  return ( x < y );
}
//-----------------------------------------------------------------------
template<typename T1, typename T2>
int ODEydot::leq( T1 x, T2 y )
{
  return ( x <= y );
}
//=======================================================================

#endif // __SBML_MATH_CPP_HPP
