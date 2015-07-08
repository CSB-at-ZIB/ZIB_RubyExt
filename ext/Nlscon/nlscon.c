#include <ruby.h>
#include "NLSCON/nlscon.h"


static VALUE gfcn;
static VALUE gjac;


void nlscon_state_free(nlscon_state_t * s)
{
  free(s->x);
  free(s->xScal);
  free(s->fi);
  free(s->fScal);
  free(s->iwk);
  free(s->rwk);

  free(s);
}

nlscon_state_t * nlscon_state_alloc(int n, int m, int mFit)
{
  int j;
  nlscon_state_t * s;

  s = (nlscon_state_t *) malloc( sizeof(nlscon_state_t) );

  s->n = n;
  s->m = m;
  s->mFit = mFit;

  s->fcn = 0;
  s->jac = 0;

      s->x = (double *) malloc( sizeof(double)*n );
  s->xScal = (double *) malloc( sizeof(double)*n );
  for (j = 0; j < n; ++j)
  {
    s->x[j] = 0.0;
    s->xScal[j] = 0.0;
  } 

     s->fi = (double *) malloc( sizeof(double)*mFit );
  s->fScal = (double *) malloc( sizeof(double)*mFit );
  for (j = 0; j < mFit; ++j)
  {
    s->fi[j] = 0.0;
    s->fScal[j] = 0.0;
  } 


  s->rTol = 1.0e-4;

  for (j = 0; j < 50; ++j)
  {
    s->iOpt[j] = 0;
  }
  s->iOpt[30] = 3;  /* nonlin: 3 */

  s->iErr = 0;

  s->lIwk = n + 52;
   s->iwk = (int *) malloc( sizeof(int)*(s->lIwk) );
  for (j = 0; j < s->lIwk; ++j)
  {
    s->iwk[j] = 0;
  }

  s->lRwk  = 2*(m + n)*n + 8*m + 10*n + ( (m > n) ? m : n ) + 104;
  s->lRwk += 2*n + n*n; /* ??? seems to be crucial to have just more memory... ??? */
   s->rwk  = (double *) malloc( sizeof(double)*(s->lRwk) );
  for (j = 0; j < s->lRwk; ++j)
  {
    s->rwk[j] = 0.0;
  }

  return s;
}

static 
void nlscon_free(void * p)
{
  nlscon_state_free(p);
}


static 
void fcn_block(int *n, int *m, int *mcon,
               double *x, double *f, int *ifail)
{
  *ifail = -987;

  {
    int j;
    VALUE arr;
    VALUE result;

    arr = rb_ary_new2(*n);

    for (j = 0; j < *n; ++j)
    {
      rb_ary_push( arr, rb_float_new(x[j]) );
    }

    result = rb_yield(arr);

    if ( TYPE(result) == T_ARRAY )
    {
      for (j = 0; j < *m; ++j)
      { 
        f[j] = NUM2DBL( rb_ary_entry(result, (long)j) );
      }

      *ifail = 0;
    }
  }
}

static
void fcn_cb(int *n, int *m, int *mcon,
            double *x, double *f, int *ifail)
{
  *ifail = -987;
 
  {
    int j;
    VALUE arr;
    VALUE result;

    arr = rb_ary_new2(*n);

    for (j = 0; j < *n; ++j)
    {
      rb_ary_push( arr, rb_float_new(x[j]) );
    }

    if ( rb_class_of(gfcn) == rb_cSymbol )
    {
      result = rb_funcall(rb_class_of(gfcn), rb_to_id(gfcn), 
                          4, INT2NUM(*n), INT2NUM(*m), INT2NUM(*mcon), arr);
    }
    else
    {
      result = rb_funcall(gfcn, rb_intern("call"), 
                          4, INT2NUM(*n), INT2NUM(*m), INT2NUM(*mcon), arr);
    }

    if ( TYPE(result) == T_ARRAY )
    {
      for (j = 0; j < *m; ++j)
      { 
        f[j] = NUM2DBL( rb_ary_entry(result, (long)j) );
      }
  
      *ifail = 0;
    }
  }
}

static 
void jac_cb(int *n, int *m, int *mcon,
            double *x, double *dfdx, int *ifail)
{
  *ifail = -987;

  {
    int j,k;
    VALUE arr;
    VALUE result;

    arr = rb_ary_new2(*n);

    for (j = 0; j < *n; ++j)
    {
      rb_ary_push( arr, rb_float_new(x[j]) );
    }

    if ( rb_class_of(gjac) == rb_cSymbol )
    {
      result = rb_funcall(rb_class_of(gjac), rb_to_id(gjac), 
                          4, INT2NUM(*n), INT2NUM(*m), INT2NUM(*mcon), arr);
    }
    else
    {
      result = rb_funcall(gjac, rb_intern("call"), 
                          4, INT2NUM(*n), INT2NUM(*m), INT2NUM(*mcon), arr);
    }


    if ( TYPE(result) == T_ARRAY )
    {
      for (j = 0; j < *m; ++j)
      { 
         arr = rb_ary_entry(result, (long)j);

         for (k = 0; k < *n; ++k)
         {
            dfdx[*m*k + j] = NUM2DBL( rb_ary_entry(arr, (long)k) );
  /* WRONG: dfdx[*n*j + k] = NUM2DBL( rb_ary_entry(arr, (long)k) ); */
         }
      }
  
      *ifail = 0;
    }
  }
}
 

static
VALUE nlscon_iterate(VALUE self)
{
  nlscon_state_t * s;

  Data_Get_Struct(self, nlscon_state_t, s);
 
  if ( rb_block_given_p() )
  {
    s->fcn = fcn_block;
  }
  
  if ( s->fcn == 0 )
  {
     rb_raise(rb_eTypeError, "Expected some callback (symbol or block).");
  }
  if ( s->jac == 0 && s->iOpt[2] == 1 )
  {
     s->iOpt[2] = 3; /* switch back to numerical differentiation (w/feedb.) */
  }

  s->iOpt[1] = 1;    /* switch stepwise mode always ON */

  nlscon_( &(s->n), &(s->m), &(s->mFit),
           s->fcn, s->jac,
           s->x, s->xScal, s->fi, s->fScal,
           &(s->rTol), s->iOpt, &(s->iErr),
           &(s->lIwk), s->iwk, &(s->lRwk), s->rwk);
 
  return INT2NUM(s->iErr);

/*
  int     j, n, m, mFit;
  double  *x, *xScal, *fi, *fScal;
  double  rTol;
  int     debug, nDense;
  int     iOpt[50];
  int     iErr, lIwk, lRwk;
  double  *iwk, *rwk;

  VALUE z;
  VALUE zscal;
  VALUE gi;
  VALUE gscal;

     n  = NUM2INT( rb_iv_get(self, "@ndim") );
     m  = NUM2INT( rb_iv_get(self, "@mdim") );
  mFit  = NUM2INT( rb_iv_get(self, "@mfit") );

      z = rb_iv_get(self, "@x");
  zscal = rb_iv_get(self, "@xscal");

      x = (double *) ALLOCA_N(double, n);
  xScal = (double *) ALLOCA_N(double, n);

  for (j = 0; j < n; ++j)
  {
         x[j] = NUM2DBL( rb_ary_entry(z, (long)j) );
     xScal[j] = NUM2DBL( rb_ary_entry(zscal, (long)j) );
  }

     gi = rb_iv_get(self, "@fobs");
  gscal = rb_iv_get(self, "@fscal");

     fi = (double *) ALLOCA_N(double, mfit);
  fScal = (double *) ALLOCA_N(double, mfit);

  for (j = 0; j < mFit; ++j)
  {
        fi[j] = NUM2DBL( rb_ary_entry(gi, (long)j) );
     fScal[j] = NUM2DBL( rb_ary_entry(gscal, (long)j) );
  }

  rTol = NUM2DBL( rb_iv_get(self, "@rtol") );

  debug = NUM2INT( rb_iv_get(self, "@monitor") );

  iOpt[0]  =  debug;    // Integration monitoring: 0 no output, 1 standard, 2 additional
  iOpt[1]  =  0;        // Unit number for monitor ( == 6 if iOpt[0] > 0 )
  iOpt[2]  =  0;        // Solution output: 0 no output, 1 initial&final vaules, 2 additional
  iOpt[3]  =  0;        // Unit number for solution ( == 6 if iOpt[2] > 0 )
  iOpt[4]  =  1;        // Singular or non-singualar matrix B: 0 sing, 1 non-sing
  iOpt[18] =
  iOpt[19] =
  iOpt[20] =
  iOpt[21] =
  iOpt[22] = -1;        // Not used in LIMEX_A (relevant in LIMEX_B, sparse Jacobians)

  iOpt[23] =            // on return: Number of function evaluations
  iOpt[24] =            // on return: Number of fcn evaltions for Jacobian computation
  iOpt[25] =            // on return: Number of LU decompositions
  iOpt[26] =            // on return: Number of back-substitions
  iOpt[27] =            // on return: Number of integration steps
  iOpt[28] =  0;        // on return: Number of Jacobian evaluations

  iOpt[29] = -1;        // Not used in LIMEX_A (relevant in LIMEX_B, sparse Jacobians)

  iErr = NUM2INT( rb_iv_get(self, "@istatus") );

  // Finally, here the routine now is called

  {
    nlscon_( &n, &m, &mFit, 
             fcn, jac, 
             &x, &xScal, &fi, &fScal, 
             &rTol, iOpt, &iErr,
             &lIwk, iwk, &lRwk, rwk
          );
  }

  VALUE y;
  VALUE yscal;

      y = rb_ary_new2(n);
  yscal = rb_ary_new2(n);

  for (j = 0; j < n; ++j)
  {
    rb_ary_push( y, rb_float_new(x[j]) );
    rb_ary_push( yscal, rb_float_new(xScal[j]) );
  }

  rb_iv_set(self, "@x", y);
  rb_iv_set(self, "@xscal", yscal);
  rb_iv_set(self, "@rtol", rb_float_new(rTol) );

  VALUE ret;
  
  ret = INT2NUM(iErr);

  rb_iv_set(self, "@istatus", ret);

  return ret;
*/
}


static
VALUE nlscon_x(VALUE self)
{
  int j;
  nlscon_state_t * s;
  VALUE arr;

  Data_Get_Struct(self, nlscon_state_t, s);

  arr = rb_ary_new2(s->n);

  for (j = 0; j < s->n; ++j)
  {
    rb_ary_push( arr, rb_float_new(s->x[j]) );
  }

  return arr;
}

static
VALUE nlscon_x_eq(VALUE self, VALUE x)
{
  int j;
  nlscon_state_t * s;

  Check_Type(x, T_ARRAY);

  Data_Get_Struct(self, nlscon_state_t, s);

  if ( s->n == RARRAY_LEN(x) )
  {
    for (j = 0; j < s->n; ++j)
    {
      s->x[j] = NUM2DBL( rb_ary_entry( x, (long)j ) );
    }   
  }
  else
  {
    rb_raise(rb_eArgError, "Wrong length in x!");
  }

  return x;
}


static
VALUE nlscon_xscal(VALUE self)
{
  int j;
  nlscon_state_t * s;
  VALUE arr;

  Data_Get_Struct(self, nlscon_state_t, s);

  arr = rb_ary_new2(s->n);

  for (j = 0; j < s->n; ++j)
  {
    rb_ary_push( arr, rb_float_new(s->xScal[j]) );
  }

  return arr;
}

static
VALUE nlscon_xscal_eq(VALUE self, VALUE xscal)
{
  int j;
  nlscon_state_t * s;

  Check_Type(xscal, T_ARRAY);

  Data_Get_Struct(self, nlscon_state_t, s);

  if ( s->n == RARRAY_LEN(xscal) )
  {
    for (j = 0; j < s->n; ++j)
    {
      s->xScal[j] = NUM2DBL( rb_ary_entry( xscal, (long)j ) );
    }   
  }
  else
  {
    rb_raise(rb_eArgError, "Wrong length in xscal!");
  }

  return xscal;
}


static
VALUE nlscon_fobs(VALUE self)
{
  int j;
  nlscon_state_t * s;
  VALUE arr;

  Data_Get_Struct(self, nlscon_state_t, s);

  arr = rb_ary_new2(s->mFit);

  for (j = 0; j < s->mFit; ++j)
  {
    rb_ary_push( arr, rb_float_new(s->fi[j]) );
  }

  return arr;
}

static
VALUE nlscon_fobs_eq(VALUE self, VALUE fobs)
{
  int j;
  nlscon_state_t * s;

  Check_Type(fobs, T_ARRAY);

  Data_Get_Struct(self, nlscon_state_t, s);

  if ( s->mFit == RARRAY_LEN(fobs) )
  {
    for (j = 0; j < s->mFit; ++j)
    {
      s->fi[j] = NUM2DBL( rb_ary_entry( fobs, (long)j ) );
    }   
  }
  else
  {
    rb_raise(rb_eArgError, "Wrong length in fobs!");
  }

  return fobs;
}


static
VALUE nlscon_fscal(VALUE self)
{
  int j;
  nlscon_state_t * s;
  VALUE arr;

  Data_Get_Struct(self, nlscon_state_t, s);

  arr = rb_ary_new2(s->mFit);

  for (j = 0; j < s->mFit; ++j)
  {
    rb_ary_push( arr, rb_float_new(s->fScal[j]) );
  }

  return arr;
}

static
VALUE nlscon_fscal_eq(VALUE self, VALUE fscal)
{
  int j;
  nlscon_state_t * s;

  Check_Type(fscal, T_ARRAY);

  Data_Get_Struct(self, nlscon_state_t, s);

  if ( s->mFit == RARRAY_LEN(fscal) )
  {
    for (j = 0; j < s->mFit; ++j)
    {
      s->fScal[j] = NUM2DBL( rb_ary_entry( fscal, (long)j ) );
    }   
  }
  else
  {
    rb_raise(rb_eArgError, "Wrong length in fscal!");
  }

  return fscal;
}


static
VALUE nlscon_iopt(VALUE self)
{
  VALUE iopt;
  nlscon_state_t * s;

  Data_Get_Struct(self, nlscon_state_t, s);  

  iopt = rb_hash_new();

  rb_hash_aset( iopt, rb_str_new2("qsucc"),  INT2NUM(s->iOpt[ 0]) );
  rb_hash_aset( iopt, rb_str_new2("mode"),   INT2NUM(s->iOpt[ 1]) );
  rb_hash_aset( iopt, rb_str_new2("jacgen"), INT2NUM(s->iOpt[ 2]) );
  rb_hash_aset( iopt, rb_str_new2("iscal"),  INT2NUM(s->iOpt[ 8]) );
  rb_hash_aset( iopt, rb_str_new2("mprerr"), INT2NUM(s->iOpt[10]) );
  rb_hash_aset( iopt, rb_str_new2("mprmon"), INT2NUM(s->iOpt[12]) );
  rb_hash_aset( iopt, rb_str_new2("mprsol"), INT2NUM(s->iOpt[14]) ); 
  /* rb_hash_aset( iopt, rb_str_new2("mprtim"), INT2NUM(s->iOpt[18]) ); */
  rb_hash_aset( iopt, rb_str_new2("qstat"),  INT2NUM(s->iOpt[20]) );
  rb_hash_aset( iopt, rb_str_new2("mprsta"), INT2NUM(s->iOpt[21]) );
  rb_hash_aset( iopt, rb_str_new2("nonlin"), INT2NUM(s->iOpt[30]) );
  rb_hash_aset( iopt, rb_str_new2("qrank1"), INT2NUM(s->iOpt[31]) );
  rb_hash_aset( iopt, rb_str_new2("qnscal"), INT2NUM(s->iOpt[34]) );
  rb_hash_aset( iopt, rb_str_new2("iterm"),  INT2NUM(s->iOpt[35]) );
  rb_hash_aset( iopt, rb_str_new2("ibdamp"), INT2NUM(s->iOpt[37]) );
  rb_hash_aset( iopt, rb_str_new2("user1"),  INT2NUM(s->iOpt[45]) );

  return iopt;
}

static
VALUE nlscon_iopt_eq(VALUE self, VALUE iopt)
{
  nlscon_state_t * s;
  VALUE val;

  Check_Type(iopt, T_HASH);

  Data_Get_Struct(self, nlscon_state_t, s);

  if ( (val=rb_hash_aref(iopt, rb_str_new2("qsucc"))) != Qnil )
  {
    s->iOpt[0] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("jacgen"))) != Qnil )
  {
    s->iOpt[2] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("iscal"))) != Qnil )
  {
    s->iOpt[8] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("mprerr"))) != Qnil )
  {
    s->iOpt[10] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("mprmon"))) != Qnil )
  {
    s->iOpt[12] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("mprsol"))) != Qnil )
  {
    s->iOpt[14] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("qstat"))) != Qnil )
  {
    s->iOpt[20] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("mprsta"))) != Qnil )
  {
    s->iOpt[21] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("nonlin"))) != Qnil )
  {
    s->iOpt[30] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("qrank1"))) != Qnil )
  {
    s->iOpt[31] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("qnscal"))) != Qnil )
  {
    s->iOpt[34] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("iterm"))) != Qnil )
  {
    s->iOpt[35] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("ibdamp"))) != Qnil )
  {
    s->iOpt[37] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iopt, rb_str_new2("user1"))) != Qnil )
  {
    s->iOpt[45] = NUM2INT( val );
  }

  return iopt;
}


static
VALUE nlscon_iwk(VALUE self)
{
  VALUE iwk;
  nlscon_state_t * s;

  Data_Get_Struct(self, nlscon_state_t, s);  

  iwk = rb_hash_new();

  rb_hash_aset( iwk, rb_str_new2("niter"),  INT2NUM(s->iwk[ 0]) );
  rb_hash_aset( iwk, rb_str_new2("ncorr"),  INT2NUM(s->iwk[ 2]) );
  rb_hash_aset( iwk, rb_str_new2("nfcn"),   INT2NUM(s->iwk[ 3]) );
  rb_hash_aset( iwk, rb_str_new2("njac"),   INT2NUM(s->iwk[ 4]) );
  rb_hash_aset( iwk, rb_str_new2("nfcnj"),  INT2NUM(s->iwk[ 7]) );
  rb_hash_aset( iwk, rb_str_new2("nrejr1"), INT2NUM(s->iwk[ 8]) );
  /* rb_hash_aset( iwk, rb_str_new2("idcode"), INT2NUM(s->iwk[11]) ); */
  rb_hash_aset( iwk, rb_str_new2("niwkfr"), INT2NUM(s->iwk[15]) );
  rb_hash_aset( iwk, rb_str_new2("nrwkfr"), INT2NUM(s->iwk[16]) );
  rb_hash_aset( iwk, rb_str_new2("liwka"),  INT2NUM(s->iwk[17]) );
  rb_hash_aset( iwk, rb_str_new2("lrwka"),  INT2NUM(s->iwk[18]) );
  rb_hash_aset( iwk, rb_str_new2("ifail"),  INT2NUM(s->iwk[22]) );
  rb_hash_aset( iwk, rb_str_new2("nitmax"), INT2NUM(s->iwk[30]) );
  rb_hash_aset( iwk, rb_str_new2("irank"),  INT2NUM(s->iwk[31]) );
  rb_hash_aset( iwk, rb_str_new2("new"),    INT2NUM(s->iwk[32]) );
  rb_hash_aset( iwk, rb_str_new2("ifccnt"), INT2NUM(s->iwk[33]) );

  return iwk;
}

static
VALUE nlscon_iwk_eq(VALUE self, VALUE iwk)
{
  nlscon_state_t * s;
  VALUE val;

  Check_Type(iwk, T_HASH);

  Data_Get_Struct(self, nlscon_state_t, s);

  if ( (val=rb_hash_aref(iwk, rb_str_new2("nitmax"))) != Qnil )
  {
    s->iwk[30] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iwk, rb_str_new2("irank"))) != Qnil )
  {
    s->iwk[31] = NUM2INT( val );
  }
  if ( (val=rb_hash_aref(iwk, rb_str_new2("ifccnt"))) != Qnil )
  {
    s->iwk[33] = NUM2INT( val );
  }

  return iwk;
}


static
VALUE nlscon_rwk(VALUE self)
{
  VALUE rwk;
  nlscon_state_t * s;

  Data_Get_Struct(self, nlscon_state_t, s);  

  rwk = rb_hash_new();

  rb_hash_aset( rwk, rb_str_new2("conv"),   rb_float_new(s->rwk[16]) );
  rb_hash_aset( rwk, rb_str_new2("sumx"),   rb_float_new(s->rwk[17]) );
  rb_hash_aset( rwk, rb_str_new2("dlevf"),  rb_float_new(s->rwk[18]) );
  rb_hash_aset( rwk, rb_str_new2("fcbnd"),  rb_float_new(s->rwk[19]) );
  rb_hash_aset( rwk, rb_str_new2("fcstrt"), rb_float_new(s->rwk[20]) );
  rb_hash_aset( rwk, rb_str_new2("fcmin"),  rb_float_new(s->rwk[21]) );
  rb_hash_aset( rwk, rb_str_new2("sigma"),  rb_float_new(s->rwk[22]) ); 
  rb_hash_aset( rwk, rb_str_new2("cond"),   rb_float_new(s->rwk[24]) );
  rb_hash_aset( rwk, rb_str_new2("ajdel"),  rb_float_new(s->rwk[25]) );
  rb_hash_aset( rwk, rb_str_new2("ajmin"),  rb_float_new(s->rwk[26]) );
  rb_hash_aset( rwk, rb_str_new2("etadif"), rb_float_new(s->rwk[27]) );
  rb_hash_aset( rwk, rb_str_new2("etaini"), rb_float_new(s->rwk[28]) );
  rb_hash_aset( rwk, rb_str_new2("prec"),   rb_float_new(s->rwk[30]) );
  rb_hash_aset( rwk, rb_str_new2("skap"),   rb_float_new(s->rwk[31]) );
  if ( s->iOpt[20] == 1 )
  {
    int j,k,n;
    VALUE xl, xr, vcv;
 
    rb_hash_aset( rwk, rb_str_new2("sigma2"), rb_float_new(s->rwk[49]) );

    n = s->n;

    xl = rb_ary_new2(n);
    xr = rb_ary_new2(n);
    vcv = rb_ary_new2(n*n);

    for (j = 0; j < n; ++j)
    {
       rb_ary_push( xl, rb_float_new(s->rwk[50+j]) );
       rb_ary_push( xr, rb_float_new(s->rwk[50+n+j]) );
       for (k = 0; k < n; ++k)
       {
          rb_ary_push( vcv, rb_float_new(s->rwk[50+2*n+n*j+k]) );
       }
    } 

    rb_hash_aset( rwk, rb_str_new2("xl"), xl );
    rb_hash_aset( rwk, rb_str_new2("xr"), xr );
    rb_hash_aset( rwk, rb_str_new2("vcv"), vcv );
  }

  return rwk;
}

static
VALUE nlscon_rwk_eq(VALUE self, VALUE rwk)
{
  nlscon_state_t * s;
  VALUE val;

  Check_Type(rwk, T_HASH);

  Data_Get_Struct(self, nlscon_state_t, s);

  if ( (val=rb_hash_aref(rwk, rb_str_new2("fcbnd"))) != Qnil )
  {
    s->rwk[19] = NUM2DBL( val );
  }
  if ( (val=rb_hash_aref(rwk, rb_str_new2("fcstrt"))) != Qnil )
  {
    s->rwk[20] = NUM2DBL( val );
  }
  if ( (val=rb_hash_aref(rwk, rb_str_new2("fcmin"))) != Qnil )
  {
    s->rwk[21] = NUM2DBL( val );
  }
  if ( (val=rb_hash_aref(rwk, rb_str_new2("sigma"))) != Qnil )
  {
    s->rwk[22] = NUM2DBL( val );
  }
  if ( (val=rb_hash_aref(rwk, rb_str_new2("cond"))) != Qnil )
  {
    s->rwk[24] = NUM2DBL( val );
  }
  if ( (val=rb_hash_aref(rwk, rb_str_new2("ajdel"))) != Qnil )
  {
    s->rwk[25] = NUM2DBL( val );
  }
  if ( (val=rb_hash_aref(rwk, rb_str_new2("ajmin"))) != Qnil )
  {
    s->rwk[26] = NUM2DBL( val );
  }
  if ( (val=rb_hash_aref(rwk, rb_str_new2("etadif"))) != Qnil )
  {
    s->rwk[27] = NUM2DBL( val );
  }
  if ( (val=rb_hash_aref(rwk, rb_str_new2("etaini"))) != Qnil )
  {
    s->rwk[28] = NUM2DBL( val );
  }

  return rwk;
}



static
VALUE nlscon_rtol_eq(VALUE self, VALUE rtol)
{
  nlscon_state_t * s;

  Check_Type(rtol, T_FLOAT);

  Data_Get_Struct(self, nlscon_state_t, s);
  s->rTol = NUM2DBL(rtol);

  return rtol;
}

static
VALUE nlscon_rtol(VALUE self)
{
  nlscon_state_t * s;

  Data_Get_Struct(self, nlscon_state_t, s);

  return rb_float_new(s->rTol);
}


static
VALUE nlscon_ierr(VALUE self)
{
  nlscon_state_t * s;

  Data_Get_Struct(self, nlscon_state_t, s);

  return INT2NUM(s->iErr);
}


static
VALUE nlscon_fcn_eq(VALUE self, VALUE fcn)
{
  nlscon_state_t * s;

  if ( (rb_class_of(fcn) != rb_cSymbol) && 
       (rb_class_of(fcn) != rb_cMethod)
     )
  {
     rb_raise(rb_eTypeError, "Expected symbol callback or method callback");
  }

  Data_Get_Struct(self, nlscon_state_t, s);

  gfcn = fcn;
  s->fcn = fcn_cb;

  /* rb_iv_set(self, "@f", fcn); */
  rb_global_variable(&gfcn);
  
  return fcn;
}

static
VALUE nlscon_jac_eq(VALUE self, VALUE jac)
{
  nlscon_state_t * s;

  if ( (rb_class_of(jac) != rb_cSymbol) &&
       (rb_class_of(jac) != rb_cMethod)
     )
  {
     rb_raise(rb_eTypeError, "Expected symbol callback or method callback");
  }

  Data_Get_Struct(self, nlscon_state_t, s);

  gjac = jac;
  s->jac = jac_cb;

  /* rb_iv_set(self, "@df", jac); */
  rb_global_variable(&gjac);
 
  return jac;
}


static 
VALUE nlscon_init(VALUE self, VALUE dims)
{
  return self;
}

 
VALUE nlscon_new(VALUE klass, VALUE dims)
{
  int n, m, mFit;
  nlscon_state_t * s;
  VALUE argv[1];
  VALUE gn;

  Check_Type(dims, T_ARRAY);

  if ( RARRAY_LEN(dims) != 3 )
  {
    rb_raise(rb_eArgError, "Wrong initial dim vector != [n, m, mFit]!");
  }

     n = NUM2INT( rb_ary_entry(dims, 0L) );
     m = NUM2INT( rb_ary_entry(dims, 1L) );
  mFit = NUM2INT( rb_ary_entry(dims, 2L) );

  s = nlscon_state_alloc(n, m, mFit);

  gn = Data_Wrap_Struct(klass, 0, nlscon_free, s);
  argv[0] = dims;

  rb_obj_call_init(gn, 1, argv);

  return gn;
}


VALUE cNlscon;

void Init_Nlscon() 
{
  cNlscon = rb_define_class("Nlscon", rb_cObject);
  rb_define_singleton_method(cNlscon, "new", nlscon_new, 1);
  rb_define_method(cNlscon, "initialize", nlscon_init, 1);
  rb_define_method(cNlscon, "iterate", nlscon_iterate, 0);
  rb_define_method(cNlscon, "x", nlscon_x, 0);
  rb_define_method(cNlscon, "x=", nlscon_x_eq, 1);
  rb_define_method(cNlscon, "xscal", nlscon_xscal, 0);
  rb_define_method(cNlscon, "xscal=", nlscon_xscal_eq, 1);
  rb_define_method(cNlscon, "fobs", nlscon_fobs, 0);
  rb_define_method(cNlscon, "fobs=", nlscon_fobs_eq, 1);
  rb_define_method(cNlscon, "fscal", nlscon_fscal, 0);
  rb_define_method(cNlscon, "fscal=", nlscon_fscal_eq, 1);
  rb_define_method(cNlscon, "rtol", nlscon_rtol, 0);
  rb_define_method(cNlscon, "rtol=", nlscon_rtol_eq, 1);
  rb_define_method(cNlscon, "iopt", nlscon_iopt, 0);
  rb_define_method(cNlscon, "iopt=", nlscon_iopt_eq, 1);
  rb_define_method(cNlscon, "ierr", nlscon_ierr, 0);
  rb_define_method(cNlscon, "iwk", nlscon_iwk, 0);
  rb_define_method(cNlscon, "iwk=", nlscon_iwk_eq, 1);
  rb_define_method(cNlscon, "rwk", nlscon_rwk, 0);
  rb_define_method(cNlscon, "rwk=", nlscon_rwk_eq, 1);
  rb_define_method(cNlscon, "f=", nlscon_fcn_eq, 1);
  rb_define_method(cNlscon, "df=", nlscon_jac_eq, 1);
}

