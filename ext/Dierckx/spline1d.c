#include <ruby.h>
#include "DIERCKX/spline1d.h"

#define MAX(a,b) ((a) > (b) ? a : b)


void spline1d_state_free(spline1d_state_t * st)
{
	free(st->t);
	free(st->coeffs);
	free(st->wrk);
	free(st->iwrk);
	
	free(st);
}

spline1d_state_t * spline1d_state_alloc(int nest, int lwrk)
{
	int j;
	spline1d_state_t * st;

	st = (spline1d_state_t *) malloc( sizeof(spline1d_state_t) );

	st->n = nest;
	st->lwrk = lwrk;

	st->t = (double*) malloc( sizeof(double)*nest );
	st->coeffs = (double*) malloc( sizeof(double)*nest );
	st->iwrk = (int*) malloc( sizeof(int)*nest );
	for (j = 0; j < nest; ++j)
	{
		st->t[j] = 0.0;
		st->coeffs[j] = 0.0;
		st->iwrk[j] = 0;
	}	

	st->wrk = (double*) malloc( sizeof(double)*lwrk );
        for (j = 0; j < lwrk; ++j)
	{
		st->wrk[j] = 0.0;
	}

	st->residual = 0.0;
	st->ierror = 0;
	
	return st;
}

static
void spline1d_free(void * p)
{
	spline1d_state_free(p);
}


static
VALUE spline1d_eval(VALUE self, VALUE x)
{
	int j, m;
	double *x_in; 
	double *y_out;
	spline1d_state_t * st;
	VALUE y;

	Check_Type(x, T_ARRAY);
	Data_Get_Struct(self, spline1d_state_t, st);

	m = RARRAY_LEN(x);
	
	x_in = (double *) malloc( sizeof(double)*m );
	y_out = (double *) malloc( sizeof(double)*m );
	for (j= 0; j < m; ++j)
	{
		x_in[j] = NUM2DBL( rb_ary_entry(x, (long)j) );
		y_out[j] = 0.0;
	}

	splev_( st->t, &(st->n), st->coeffs, &(st->k), 
		x_in, y_out, &m, &(st->bc), &(st->ierror) );

	y = rb_ary_new2(m);
	for (j = 0; j < m; ++j)
	{
		rb_ary_push( y, rb_float_new(y_out[j]) );
	}

	return y;
}

static
VALUE spline1d_derivative(VALUE self, VALUE nu, VALUE x)
{
	int j, m, nu_in;
	double *x_in; 
	double *y_out;
	spline1d_state_t * st;
	VALUE y;

	Check_Type(x, T_ARRAY);
	Data_Get_Struct(self, spline1d_state_t, st);

	m = RARRAY_LEN(x);
	nu_in = NUM2INT(nu);

	x_in = (double *) malloc( sizeof(double)*m );
	y_out = (double *) malloc( sizeof(double)*m );
	for (j= 0; j < m; ++j)
	{
		x_in[j] = NUM2DBL( rb_ary_entry(x, (long)j) );
		y_out[j] = 0.0;
	}

	splder_( st->t, &(st->n), st->coeffs, &(st->k), &nu_in,
		x_in, y_out, &m, &(st->bc), st->wrk, &(st->ierror) );

	y = rb_ary_new2(m);
	for (j = 0; j < m; ++j)
	{
		rb_ary_push( y, rb_float_new(y_out[j]) );
	}

	return y;
}


static
VALUE spline1d_show(VALUE self)
{
  VALUE info;
  spline1d_state_t * st;

  Data_Get_Struct(self, spline1d_state_t, st);  

  info = rb_hash_new();

  rb_hash_aset( info, rb_str_new2("n"),         INT2NUM(st->n) );
  rb_hash_aset( info, rb_str_new2("k"),         INT2NUM(st->k) );
  rb_hash_aset( info, rb_str_new2("bc"),        INT2NUM(st->bc) );
  rb_hash_aset( info, rb_str_new2("bc_text"),   rb_str_new2(
			  				(st->bc==0) ? "extrapolate" :
			 				(st->bc==1) ? "zero" :
			 				(st->bc==2) ? "error" :
			 				(st->bc==3) ? "nearest" : "unknown bc") );
  rb_hash_aset( info, rb_str_new2("residual"),  rb_float_new(st->residual) );

  return info;
}


static
VALUE spline1d_ierror(VALUE self)
{
  spline1d_state_t * st;

  Data_Get_Struct(self, spline1d_state_t, st);

  return INT2NUM(st->ierror);
}



static
VALUE spline1d_init(int argc, VALUE*argv, VALUE self)
{
	return self;
}


VALUE spline1d_new(int argc, VALUE *argv, VALUE klass)
{
  int j, m_in;
  double *x_in, *y_in, *w_in;
  double s_in;
  int k_in, bc_in, periodic;
  VALUE x, y;
  VALUE k, s, bc, period;
  VALUE spl;

  rb_scan_args( argc, argv, "24", 
		&x, &y, 
		&k, &s, &bc, &period );


  periodic = (!NIL_P(period))? RTEST(period) : 0;
  bc_in = (!NIL_P(bc)) ? NUM2INT(bc) : 3;
  s_in =  (!NIL_P(s)) ? NUM2DBL(s) : 0.0;
  k_in =  (!NIL_P(k)) ? NUM2INT(k) : 3;


  Check_Type(x, T_ARRAY);
  Check_Type(y, T_ARRAY);

  m_in = RARRAY_LEN(x);

  if ( RARRAY_LEN(x) != RARRAY_LEN(y) )
  {
    rb_raise(rb_eArgError, "Mismatch of vector sizes: len(x) != len(y)");
  }
  if ( m_in <= k_in )
  {
    rb_raise(rb_eArgError, "Order k must be less than len(x)");
  }
  if (k_in < 1 || k_in > 5)
  {
    rb_raise(rb_eArgError, "Order k must hold 1 <= k <= 5");
  }

  x_in = (double *) malloc( sizeof(double)*m_in );
  y_in = (double *) malloc( sizeof(double)*m_in );
  w_in = (double *) malloc( sizeof(double)*m_in );

  for (j= 0; j < m_in; ++j)
  {
	  x_in[j] = NUM2DBL( rb_ary_entry(x, (long)j ));
	  y_in[j] = NUM2DBL( rb_ary_entry(y, (long)j ));
	  w_in[j] = 1.0;
  }


  int iopt = 0;
  int nest = (periodic) 
	     ? MAX(m_in + 2*k_in, 2*k_in + 3) 
	     : MAX(m_in +k_in + 1, 2*k_in + 3);
  int lwrk = (periodic)
	     ? m_in * (k_in+1) + nest*(8 + 5*k_in)
	     : m_in * (k_in+1) + nest*(7 + 3*k_in);

  spline1d_state_t * st = spline1d_state_alloc(nest, lwrk);

  st->bc = bc_in;
  st->k = k_in;

  if (!periodic)
  {
	  double xb = x_in[0];
	  double xe = x_in[m_in-1];

	  curfit_( &iopt, &m_in, x_in, y_in, w_in, 
		   &xb, &xe, &k_in, &s_in, &nest, 
		   &(st->n), st->t, st->coeffs, &(st->residual),
		   st->wrk, &(st->lwrk), st->iwrk, &(st->ierror) );
  }
  else
  {
	  percur_( &iopt, &m_in, x_in, y_in, w_in,
		   &k_in, &s_in, &nest,
		   &(st->n), st->t, st->coeffs, &(st->residual),
		   st->wrk, &(st->lwrk), st->iwrk, &(st->ierror) );
  }



  spl = Data_Wrap_Struct(klass, 0, spline1d_free, st);

  rb_obj_call_init(spl, argc, argv);

  return spl;
}


VALUE cSpline1D;

void Init_Dierckx() 
{
  cSpline1D = rb_define_class("Spline1D", rb_cObject);
  rb_define_singleton_method(cSpline1D, "new", spline1d_new, -1);
  rb_define_method(cSpline1D, "initialize", spline1d_init, -1);
  rb_define_method(cSpline1D, "eval", spline1d_eval, 1);
  rb_define_method(cSpline1D, "derivative", spline1d_derivative, 2);
  /*
  rb_define_method(cSpline1D, "integrate", spline1d_integrate, 3);
  rb_define_method(cSpline1D, "roots", spline1d_roots, 1);
  rb_define_method(cSpline1D, "knots", spline1d_knots, 0);
  rb_define_method(cSpline1D, "coeffs", spline1d_coeffs, 0);
  rb_define_method(cSpline1D, "residual", spline1d_residual, 0);
  */
  rb_define_method(cSpline1D, "ierror", spline1d_ierror, 0);
  rb_define_method(cSpline1D, "show", spline1d_show, 0);
}


