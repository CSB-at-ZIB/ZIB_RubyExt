#include <ruby.h>
#include "LIMEX4_3A/LIMEX4_3A.h"


VALUE cLimex;

/*
static 
void limex_free(void *p)
{
}
*/

static
void eval_hermite(int n, double t, double *y, 
                  int kOrder, double *coeff, double t1, double t2) 
{
  int j,k;
  double tFrac, tmp; 

  if (t1 == t2) return;

  tFrac = (t - t1) / (t2 - t1);
  tmp = tFrac - 1.0;

  for (j = 0; j < n; ++j)
  {
    y[j] = coeff[MAX_NO_EQNS*(kOrder+1) + j];
  }

  for (k = 0; k < kOrder; ++k)
  {
    for (j = 0; j < n; ++j)
    {
      y[j] *= tmp;
      y[j] += coeff[MAX_NO_EQNS*(kOrder-k) + j];
    }
  }

  for (j = 0; j < n; ++j)
  {
    y[j] *= tFrac;
    y[j] += coeff[j];
  }
}


static 
void fcn(int *n, int *nz,
         double *t, double *y, double *dy,
         double *B, int *ir, int *ic,
         int *info)
{
  *info = -987;

  if ( rb_block_given_p() )
  {
    int j;
    VALUE arr, z;
    VALUE result;

    z = rb_ary_new2(*n);

    for (j = 0; j < *n; ++j)
    {
      rb_ary_push( z, rb_float_new(y[j]) );
    }

    arr = rb_ary_new();
    rb_ary_push( arr, rb_float_new(*t) );
    rb_ary_push( arr, z );
    /* rb_ary_push( arr, rb_ary_new2(*n) ); */

    result = rb_yield(arr);

    if ( TYPE(result) == T_ARRAY )
    {
      for (j = 0; j < *n; ++j)
      { 
        dy[j] = NUM2DBL( rb_ary_entry(result, (long)j) );
      }

      *nz = *n;

      for (j = 0; j < *nz; ++j)
      {
          B[j] = 1.0;
         ir[j] = ic[j] = j+1;
      }

      *info = 0;
    }
  }
}


static
VALUE limex_srun(VALUE self, VALUE tspan, VALUE pidx)
{
  Check_Type(tspan, T_ARRAY);
  Check_Type(pidx, T_ARRAY);

  int     j, k, n, nDAE, nTp, nPdx;
  double  t0, T;
  double  *tp, *z, *dz;
  double  rtol, atol;
  double  h, hMax, hInit;
  int     debug, nDense;
  int     iOpt[32];
  double  rOpt[5];
  int     *iPos, iFail[3];
  int     kOrder;
  double  dense[MAX_NO_EQNS*(2 + MAX_ROW_TAB*(MAX_ROW_TAB+1)/2)];
  double  t1, t2;
  double  *y, *rTol, *aTol;

  VALUE y0;
  VALUE steps;
  VALUE solution;

  nTp = RARRAY_LEN(tspan);
  nPdx = RARRAY_LEN(pidx);

  if ( nTp < 2 )
  {
    rb_raise(rb_eArgError, "tspan array length wrong (must be 2, at least).");
    return Qnil;
  }
  if ( nPdx < 1 )
  {
    rb_raise(rb_eArgError, "pidx array length wrong (must be greater than 0).");
    return Qnil;
  }

  // rb_iv_set(self, "@t0", tStart);
  nDAE = NUM2INT( rb_iv_get(self, "@dim") );
  n    = nDAE*(1 + nPdx);
  // t0 = NUM2DBL(tStart);
  // T  = NUM2DBL(tEnd);

  rtol = NUM2DBL( rb_iv_get(self, "@rtol") );
  atol = NUM2DBL( rb_iv_get(self, "@atol") );
    y0 = rb_iv_get(self, "@y0");

    tp = (double *) ALLOCA_N(double, nTp);
    y  = (double *) ALLOCA_N(double, n);
    z  = (double *) ALLOCA_N(double, n);
    dz = (double *) ALLOCA_N(double, n);
  rTol = (double *) ALLOCA_N(double, n);
  aTol = (double *) ALLOCA_N(double, n);

  for (k = 0; k < nTp; ++k)
  {
    tp[k] = NUM2DBL( rb_ary_entry(tspan, (long)k) );
  }
  t0 = tp[0];
  T = tp[nTp-1];

  for (j = 0; j < n; ++j)
  {
     if ( rtol < 1.0e-6 )  // 1.0e-5 ?!
     {
        rTol[j] = (j < nDAE) ? rtol : 100.0*rtol;
     } 
     else
     {
        rTol[j] = (j < nDAE) ? rtol : 1.0e-4;
     }

     if ( atol < 1.0e-6 )  // 1.0e-5 ?!
     {
        aTol[j] = (j < nDAE) ? atol : 100.0*atol;
     }
     else
     {
        aTol[j] = (j < nDAE) ? atol : 1.0e-4;
     }

      z[j] = (j < nDAE) ? NUM2DBL( rb_ary_entry(y0, (long)j) ) : 0.0;
     dz[j] = 0.0;
  }

  hMax = NUM2DBL( rb_iv_get(self, "@hmax") );
  hInit = NUM2DBL( rb_iv_get(self, "@inistep") );
  h = hInit /* rTol */;

  debug = NUM2INT( rb_iv_get(self, "@monitor") );
  nDense = NUM2INT( rb_iv_get(self, "@dense") );

  iOpt[0]  =  debug;    // Integration monitoring: 0 no output, 1 standard, 2 additional
  iOpt[1]  =  0;        // Unit number for monitor ( == 6 if iOpt[0] > 0 )
  iOpt[2]  =  0;        // Solution output: 0 no output, 1 initial&final vaules, 2 additional
  iOpt[3]  =  0;        // Unit number for solution ( == 6 if iOpt[2] > 0 )
  iOpt[4]  =  1;        // Singular or non-singualar matrix B: 0 sing, 1 non-sing
  iOpt[5]  =  0;        // Determination of consistent initial values (CIV): 0 no, 1 determ
  iOpt[6]  =  0;        // Numerical or analytical Jacobian: 0 num diff approx, 1 analytical

  iOpt[7]  =  n;        // Lower bandwidth of Jacobian: 0 <= iOpt[7] <= n <= Max_Lower_Diags
  iOpt[8]  =  n;        // Upper bandwidth of Jacobian: 0 <= iOpt[8] <= n <= Max_Upper_Diags

  iOpt[9]  =  1;        // Re-use of Jacobian: 0 no re-use, 1 re-use of Jacobian in the following steps
  iOpt[10] =  1;        // Switch for error tolerances: 0 rTol&aTol scalar, 1 rTol&aTol are vectors
  iOpt[11] =  1;        // Switch for one step mode: 0 off, 1 return from each step, 2 return only from prescribed steps

  iOpt[12] =  0;        // Dense output option: 0 off, 1 on equidist pts within interval, 2 on equidist pts within step (# in iOpt[13]), 3 on additional pts
  iOpt[13] =  0;        // Number of equidistant points if iOpt[12] == 1 or 2
  iOpt[14] =  0;        // Unit number for dense output (iOpt[14] == 0 suppresses dense output)

  iOpt[15] =  0;        // Type of call, may be modified! 0 initial call, 1 successive call

  iOpt[16] =  0;        // Behaviour at t_End: 0 stop exactly at t_End, 1 may compute/use values also for t > t_End
  iOpt[17] =  0;        // PostScript plot of Jacobian: 0 no plot, j plot at j-th step, -1 plot for initial step (step 0)

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

  iOpt[30] = -1;        // !!! Only available in LIMD !!!
                        // Type of left-hand side B: 0 B=id, 1 B=const., 2 variable B
  iOpt[31] = -1;        // !!! Only available in (s)LIMDHERM !!!
                        // Interpolation mode: 0 no additional output, 1 give additional output (switched on)

  if ( (nTp == 2) && (nDense > 0) )
  {
    iOpt[11] = 2;

    iOpt[12] = 1;
    iOpt[13] = nDense;
  }

  rOpt[0] = rOpt[1] = rOpt[2] /* = rOpt[3] = rOpt[4] */ = 0.0;
  rOpt[0] = hMax;

  iPos = (int *) ALLOCA_N(int, n);

  for (j = 0; j < n; ++j)
  {
    iPos[j] = 0;
  }

  iFail[0] = iFail[1] = iFail[2] = 0;

  VALUE tmparr;

  tmparr = rb_ary_new2(n);
  steps = rb_ary_new();
  solution = rb_hash_new();


  for (j = 0; j < n; ++j)
  {
    rb_ary_push( tmparr, rb_float_new(z[j]) ); // z is loaded with (y0,0...0)
  }

  rb_ary_push( steps, rb_float_new(t0) );
  rb_hash_aset( solution, rb_float_new(t0), tmparr ); 


  // Finally, here starts the central 
  //      *single step w/hermite* sLIMDHERM loop (iOpt[31]==1)
  //  or  *single step*           sLIMDHERM loop (iOpt[11]==1)
  //  or  *dense output*          sLIMDHERM loop (iOpt[11]==2, iOpt[13]==nDense)

  if ( nTp > 2 )
  {
    iOpt[31] = 1;  // switch on interpolation mode in sLIMDHERM

    k = 1;

      while ( (iFail[0] == 0) && (t0 < T) )
      {
        slimdherm_( &nDAE, &n, fcn, 0, &t0, &T, z, dz, 
                    rTol, aTol, &h, iOpt, rOpt, iPos, 
                    iFail, &kOrder, dense, &t1, &t2
                  );

        // fprintf(stderr,"k=%d, t0=%f, t1=%f, t2=%f, tp[k]=%f\n",k,t0,t1,t2,tp[k]);

        if (t0 <= T)
        {
          if (t1 == t2)
          {
            while ( (k < nTp) && (tp[k] <= t2) )
            {
              VALUE stp;
              VALUE arr;

              stp = rb_float_new(tp[k]);
              arr = rb_ary_new2(n);

              rb_ary_push( steps, stp );

              for (j = 0; j < n; ++j)
              {
                rb_ary_push( arr, rb_float_new(z[j]) );
              }

              rb_hash_aset( solution, stp, arr );

              ++k;
            }
          }
          else
          {
            while ( (k < nTp) && (tp[k] <= t2) )
            {
              VALUE stp;
              VALUE arr;

              stp = rb_float_new(tp[k]);
              arr = rb_ary_new2(n);

              rb_ary_push( steps, stp );

              eval_hermite(n, tp[k], y, kOrder, dense, t1, t2);

              for (j = 0; j < n; ++j)
              {
                rb_ary_push( arr, rb_float_new(y[j]) );
              }

              rb_hash_aset( solution, stp, arr );

              ++k;
            }
          }
        }
      }
  }
  else
  {
    while ( (iFail[0] == 0) && (t0 < T) )
    {
      slimdherm_( &nDAE, &n, fcn, 0, &t0, &T, z, dz, 
                  rTol, aTol, &h, iOpt, rOpt, iPos, 
                  iFail, &kOrder, dense, &t1, &t2
                );

      if (t0 <= T)
      {
        VALUE stp;
        VALUE arr;

        stp = rb_float_new(t0);
        arr = rb_ary_new2(n);

        rb_ary_push( steps, stp );

        for (j = 0; j < n; ++j)
        {
          rb_ary_push( arr, rb_float_new(z[j]) );
        }

        rb_hash_aset( solution, stp, arr );
      }
    }
  }

  rb_iv_set(self, "@t0", rb_float_new(tp[0]));
  rb_iv_set(self, "@T", rb_float_new(tp[nTp-1]));
  rb_iv_set(self, "@steps", steps);
  rb_iv_set(self, "@solution", solution);

  VALUE arr;
  arr = rb_ary_new();
  rb_ary_push( arr, INT2NUM(iFail[0]) );
  rb_ary_push( arr, INT2NUM(iFail[1]) );
  rb_ary_push( arr, INT2NUM(iFail[2]) );

  rb_iv_set(self, "@ifail", arr);

  return arr;
}


static
VALUE limex_run(VALUE self, VALUE tspan)
{
  Check_Type(tspan, T_ARRAY);

  int     j, k, n, nTp;
  double  t0, T;
  double  *tp, *z, *dz;
  double  rTol, aTol;
  double  h, hMax, hInit;
  int     debug, nDense;
  int     iOpt[32];
  double  rOpt[5];
  int     *iPos, iFail[3];
  int     kOrder;
  double  dense[MAX_NO_EQNS*(2 + MAX_ROW_TAB*(MAX_ROW_TAB+1)/2)];
  double  t1, t2;
  double  *y;

  VALUE y0;
  VALUE steps;
  VALUE solution;

  nTp = RARRAY_LEN(tspan);

  if ( nTp < 2 )
  {
    rb_raise(rb_eArgError, "tspan array length wrong (must be 2, at least).");
    return Qnil;
  }

  // rb_iv_set(self, "@t0", tStart);
  n  = NUM2INT( rb_iv_get(self, "@dim") );
  // t0 = NUM2DBL(tStart);
  // T  = NUM2DBL(tEnd);

  y0 = rb_iv_get(self, "@y0");

  tp = (double *) ALLOCA_N(double, nTp);
  y  = (double *) ALLOCA_N(double, n);
  z  = (double *) ALLOCA_N(double, n);
  dz = (double *) ALLOCA_N(double, n);

  for (k = 0; k < nTp; ++k)
  {
    tp[k] = NUM2DBL( rb_ary_entry(tspan, (long)k) );
  }
  t0 = tp[0];
  T = tp[nTp-1];

  for (j = 0; j < n; ++j)
  {
     z[j] = NUM2DBL( rb_ary_entry(y0, (long)j) );
    dz[j] = 0.0;
  }

  rTol = NUM2DBL( rb_iv_get(self, "@rtol") );
  aTol = NUM2DBL( rb_iv_get(self, "@atol") );
  hMax = NUM2DBL( rb_iv_get(self, "@hmax") );
  hInit = NUM2DBL( rb_iv_get(self, "@inistep") );
  h = hInit /* rTol */;

  debug = NUM2INT( rb_iv_get(self, "@monitor") );
  nDense = NUM2INT( rb_iv_get(self, "@dense") );

  iOpt[0]  =  debug;    // Integration monitoring: 0 no output, 1 standard, 2 additional
  iOpt[1]  =  0;        // Unit number for monitor ( == 6 if iOpt[0] > 0 )
  iOpt[2]  =  0;        // Solution output: 0 no output, 1 initial&final vaules, 2 additional
  iOpt[3]  =  0;        // Unit number for solution ( == 6 if iOpt[2] > 0 )
  iOpt[4]  =  1;        // Singular or non-singualar matrix B: 0 sing, 1 non-sing
  iOpt[5]  =  0;        // Determination of consistent initial values (CIV): 0 no, 1 determ
  iOpt[6]  =  0;        // Numerical or analytical Jacobian: 0 num diff approx, 1 analytical

  iOpt[7]  =  n;        // Lower bandwidth of Jacobian: 0 <= iOpt[7] <= n <= Max_Lower_Diags
  iOpt[8]  =  n;        // Upper bandwidth of Jacobian: 0 <= iOpt[8] <= n <= Max_Upper_Diags

  iOpt[9]  =  1;        // Re-use of Jacobian: 0 no re-use, 1 re-use of Jacobian in the following steps
  iOpt[10] =  0;        // Switch for error tolerances: 0 rTol&aTol scalar, 1 rTol&aTol are vectors
  iOpt[11] =  1;        // Switch for one step mode: 0 off, 1 return from each step, 2 return only from prescribed steps

  iOpt[12] =  0;        // Dense output option: 0 off, 1 on equidist pts within interval, 2 on equidist pts within step (# in iOpt[13]), 3 on additional pts
  iOpt[13] =  0;        // Number of equidistant points if iOpt[12] == 1 or 2
  iOpt[14] =  0;        // Unit number for dense output (iOpt[14] == 0 suppresses dense output)

  iOpt[15] =  0;        // Type of call, may be modified! 0 initial call, 1 successive call

  iOpt[16] =  0;        // Behaviour at t_End: 0 stop exactly at t_End, 1 may compute/use values also for t > t_End
  iOpt[17] =  0;        // PostScript plot of Jacobian: 0 no plot, j plot at j-th step, -1 plot for initial step (step 0)

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

  iOpt[30] = -1;        // !!! Only available in LIMD !!!
                        // Type of left-hand side B: 0 B=id, 1 B=const., 2 variable B
  iOpt[31] = -1;        // !!! Only available in LIMDHERM !!!
                        // Interpolation mode: 0 no additional output, 1 give additional output (switched on)

  if ( (nTp == 2) && (nDense > 0) )
  {
    iOpt[11] = 2;

    iOpt[12] = 1;
    iOpt[13] = nDense;
  }

  rOpt[0] = rOpt[1] = rOpt[2] /* = rOpt[3] = rOpt[4] */ = 0.0;
  rOpt[0] = hMax;

  iPos = (int *) ALLOCA_N(int, n);

  for (j = 0; j < n; ++j)
  {
    iPos[j] = 0;
  }

  iFail[0] = iFail[1] = iFail[2] = 0;

  steps = rb_ary_new();
  solution = rb_hash_new();

  rb_ary_push( steps, rb_float_new(t0) );
  rb_hash_aset( solution, rb_float_new(t0), y0 );

  // Finally, here starts the central 
  //      *single step w/hermite* LIMDHERM loop (iOpt[31]==1)
  //  or  *single step*           LIMEX    loop (iOpt[11]==1)
  //  or  *dense output*          LIMEX    loop (iOpt[11]==2, iOpt[13]==nDense)

  if ( nTp > 2 )
  {
    iOpt[31] = 1;  // switch on interpolation mode in LIMDHERM

    k = 1;
    //!!! old integration along given tp array !!!
    //for (k = 1; k < nTp; ++k)
    //{
    //  t0 = tp[k-1];
    //  T = tp[k];

      while ( (iFail[0] == 0) && (t0 < T) )
      {
        limdherm_( &n, fcn, 0, &t0, &T, z, dz, 
                   &rTol, &aTol, &h, iOpt, rOpt, iPos, 
                   iFail, &kOrder, dense, &t1, &t2
              );

        // fprintf(stderr,"k=%d, t0=%f, t1=%f, t2=%f, tp[k]=%f\n",k,t0,t1,t2,tp[k]);

        if (t0 <= T)
        {
          if (t1 == t2)
          {
            while ( (k < nTp) && (tp[k] <= t2) )
            {
              VALUE stp;
              VALUE arr;

              stp = rb_float_new(tp[k]);
              arr = rb_ary_new2(n);

              rb_ary_push( steps, stp );

              for (j = 0; j < n; ++j)
              {
                rb_ary_push( arr, rb_float_new(z[j]) );
              }

              rb_hash_aset( solution, stp, arr );

              ++k;
            }
          }
          else
          {
            while ( (k < nTp) && (tp[k] <= t2) )
            {
              VALUE stp;
              VALUE arr;

              stp = rb_float_new(tp[k]);
              arr = rb_ary_new2(n);

              rb_ary_push( steps, stp );

              eval_hermite(n, tp[k], y, kOrder, dense, t1, t2);

              for (j = 0; j < n; ++j)
              {
                rb_ary_push( arr, rb_float_new(y[j]) );
              }

              rb_hash_aset( solution, stp, arr );

              ++k;
            }
          }
        }

        /*
        if ( (iFail[0] == 0) && (t0 == T) )
        {
          VALUE stp;
          VALUE arr;

          stp = rb_float_new(t0);
          arr = rb_ary_new2(n);

          rb_ary_push( steps, stp );

          for (j = 0; j < n; ++j)
          {
            rb_ary_push( arr, rb_float_new(z[j]) );
          }

          rb_hash_aset( solution, stp, arr );
        }

        if ( iFail[0] != 0 )
        {
          VALUE stp;
          VALUE arr;

          stp = rb_float_new(t0);
          arr = rb_ary_new2(n);

          rb_ary_push( steps, stp );

          for (j = 0; j < n; ++j)
          {
            rb_ary_push( arr, rb_float_new(0.0) );
          }

          rb_hash_aset( solution, stp, arr );
        }
        */
      }
    //!!!old for loop over tp[k]
    //}
  }
  else
  {
    while ( (iFail[0] == 0) && (t0 < T) )
    {
      limex_( &n, fcn, 0, &t0, &T, z, dz, 
              &rTol, &aTol, &h, iOpt, rOpt, iPos, 
              iFail
            );

      if (t0 <= T)
      {
        VALUE stp;
        VALUE arr;

        stp = rb_float_new(t0);
        arr = rb_ary_new2(n);

        rb_ary_push( steps, stp );

        for (j = 0; j < n; ++j)
        {
          rb_ary_push( arr, rb_float_new(z[j]) );
        }

        rb_hash_aset( solution, stp, arr );
      }
    }
  }

  rb_iv_set(self, "@t0", rb_float_new(tp[0]));
  rb_iv_set(self, "@T", rb_float_new(tp[nTp-1]));
  rb_iv_set(self, "@steps", steps);
  rb_iv_set(self, "@solution", solution);

  VALUE arr;
  arr = rb_ary_new();
  rb_ary_push( arr, INT2NUM(iFail[0]) );
  rb_ary_push( arr, INT2NUM(iFail[1]) );
  rb_ary_push( arr, INT2NUM(iFail[2]) );

  rb_iv_set(self, "@ifail", arr);

  return arr;
}


static
VALUE limex_interval(VALUE self)
{
  VALUE arr;

  arr = rb_ary_new();
  rb_ary_push( arr, rb_iv_get(self, "@t0") );
  rb_ary_push( arr, rb_iv_get(self, "@T") );

  return arr;
}


static
VALUE limex_steps(VALUE self)
{
  return rb_iv_get(self, "@steps");
}


static
VALUE limex_solution(VALUE self)
{
  return rb_iv_get(self, "@solution");
}


static
VALUE limex_ifail(VALUE self)
{
  return rb_iv_get(self, "@ifail");
}

static
VALUE limex_status(VALUE self)
{
  VALUE ifail;
  char  *msg = "n/a";

  ifail = rb_iv_get(self, "@ifail");

  if ( RARRAY_LEN(ifail) == 3 )
  {
    int rc = NUM2INT( rb_ary_entry(ifail, 0L) );

    switch (rc)
    {
       case   0 : msg = "ok";
                  break;
       case  -1 : msg = "ODE dim out of range (Max_Nr_of_Equations)";
                  break;
       case  -2 : msg = "rTol(j) <= 0 for some index j";
                  break;
       case  -3 : msg = "aTol(j) < 0 for some index j";
                  break;
       case  -4 : msg = "internal range error:  Iopt(1) < 0  or  Iopt(1) > 2  (integration monitor)";
                  break;
       case  -5 : msg = "internal range error:  Iopt(2) < 0  and  Iopt(1) > 0  (monitor output unit)";
                  break;
       case  -6 : msg = "internal range error:  Iopt(3) < 0  or  Iopt(3) > 2  (solution output)";
                  break;
       case  -7 : msg = "internal range error:  Iopt(4) < 0  and  Iopt(3) > 0  (solution output unit)";
                  break;
       case  -8 : msg = "internal range error:  Iopt(5) < 0  or  Iopt(5) > 1  (matrix B singular/nonsing.)";
                  break;
       case  -9 : msg = "internal range error:  Iopt(6) < 0  or  Iopt(6) > 1  (CIV compuation)";
                  break;
       case -10 : msg = "internal range error:  Iopt(7) < 0  or  Iopt(7) > 1  (numeric/analytic Jacobian)";
                  break;
       case -11 : msg = "internal range error:  Iopt(8) out of range (Max_Lower_Diagonals)";
                  break;
       case -12 : msg = "internal range error:  Iopt(9) out of range (Max_Upper_Diagonals)";
                  break;
       case -13 : msg = "internal range error:  Iopt(10) < 0  or  Iopt(10) > 1  (reuse of Jacobian)";
                  break;
       case -14 : msg = "internal range error:  Iopt(11) < 0  or  Iopt(11) > 1  (scalar/vector rTol,aTol)";
                  break;
       case -15 : msg = "internal range error:  Iopt(12) < 0  or  Iopt(12) > 1  (one step mode)";
                  break;
       case -16 : msg = "internal range error:  Iopt(13) < 0  or  Iopt(13) > 3  (dense output)";
                  break;
       case -17 : msg = "internal range error:  Iopt(14) < 0  and  Iopt(13) = 1 or 2  (number equidist. points)";
                  break;
       case -18 : msg = "internal range error:  Iopt(15) < 0  and  Iopt(13) > 0  (dense output unit)";
                  break;
       case -19 : msg = "internal range error:  Iopt(16) < 0  or  Iopt(16) > 1  (type of call)";
                  break;
       case -20 : msg = "internal range error:  Iopt(17) < 0  or  Iopt(17) > 1  (integration exactly up to t_End)";
                  break;
       case -21 : msg = "internal range error:  Iopt(18) < -1  (postscript plot(s) of Jacobian)";
                  break;

       case -27 : msg = "internal range error:  Ropt(1) < 0  (maximal allowed stepsize)";
                  break;
       case -28 : msg = "internal range error:  Ropt(2) <= 0  and  Iopt(13) = 3  (maximal distance of dense output points)";
                  break;

       case -32 : msg = "an initial value y(j) < 0, but according to IPos it should be >= 0";
                  break;
       case -33 : msg = "user routine FCN returns an error in the first call of the current integration";
                  break;
       case -34 : msg = "more than Max_Non_Zeros_B non-zero entries in B(t,y) defined";
                  break;
       case -35 : msg = "user routine FCN returns an error in the numerical evaluation of the Jacobian";
                  break;
       case -36 : msg = "user routine JAC returns an error";
                  break;

       case -39 : msg = "internal error in LAPACK dgetrf or dgbtrf; see LAPACK User's Guide";
                  break;

       case -43 : msg = "problem not solvable with this LIMEX version; probably index of DAE > 1";
                  break;
       case -44 : msg = "problem not solvable with this LIMEX version; probably initial values not consistent or index of DAE > 1";
                  break;
       case -45 : msg = "during CIV computation, some solution y(j) is negative, but according to IPos it should be >= 0";
                  break;
       case -46 : msg = "more stepsize reductions than allowed (parameter: Max_Step_Red_Ex)";
                  break;
       case -47 : msg = "singular matrix pencil B - hA; problem not solvable with LIMEX";
                  break;
       case -48 : msg = "more integration steps than allowed (parameter: Max_Int_Steps)";
                  break;
       case -49 : msg = "more internal Newton steps for CIV computation than allowed (parameter: Max_Step_CIV)";
                  break;
       case -50 : msg = "stepsize too small; most probably too many stepsize reductions";
                  break;

        default : msg = "unknown error";
                  break;
    }
  }
  else
  {
    msg = "n/a (no solver calls so far)";
  } 

  return rb_str_new2(msg);
}


static
VALUE limex_y0(VALUE self)
{
  return rb_iv_get(self, "@y0");
}

static
VALUE limex_y0_eq(VALUE self, VALUE y0)
{
  Check_Type(y0, T_ARRAY);

  rb_iv_set(self,  "@dim",  INT2NUM( RARRAY_LEN(y0) ) );
  rb_iv_set(self,  "@y0",   y0);

  return y0;
}


static
VALUE limex_monitor(VALUE self)
{
  return rb_iv_get(self, "@monitor");
}

static
VALUE limex_monitor_eq(VALUE self, VALUE dbg)
{
  rb_iv_set(self,  "@monitor",   dbg);

  return dbg;
}


static
VALUE limex_tolerance(VALUE self)
{
  VALUE tols;
  
  tols = rb_hash_new();

  rb_hash_aset( tols, rb_str_new2("rtol"),    rb_iv_get(self, "@rtol") );
  rb_hash_aset( tols, rb_str_new2("atol"),    rb_iv_get(self, "@atol") );
  rb_hash_aset( tols, rb_str_new2("inistep"), rb_iv_get(self, "@inistep") );
  rb_hash_aset( tols, rb_str_new2("hmax"),    rb_iv_get(self, "@hmax") );

  return tols;
}

static
VALUE limex_rtol_eq(VALUE self, VALUE rtol)
{
  Check_Type(rtol, T_FLOAT);

  rb_iv_set(self, "@rtol", rtol);

  return rtol;
}

static
VALUE limex_atol_eq(VALUE self, VALUE atol)
{
  Check_Type(atol, T_FLOAT);

  rb_iv_set(self, "@atol", atol);

  return atol;
}


static
VALUE limex_inistep(VALUE self)
{
  return rb_iv_get(self, "@inistep");
}

static
VALUE limex_inistep_eq(VALUE self, VALUE inistep)
{
  Check_Type(inistep, T_FLOAT);

  rb_iv_set(self, "@inistep", inistep);

  return inistep;
}


static
VALUE limex_hmax(VALUE self)
{
  return rb_iv_get(self, "@hmax");
}

static
VALUE limex_hmax_eq(VALUE self, VALUE hmax)
{
  Check_Type(hmax, T_FLOAT);

  rb_iv_set(self, "@hmax", hmax);

  return hmax;
}


static
VALUE limex_dense(VALUE self)
{
  return rb_iv_get(self, "@dense");
}

static
VALUE limex_dense_eq(VALUE self, VALUE dense)
{
  return rb_iv_set(self, "@dense", dense);
}


static 
VALUE limex_init(VALUE self, VALUE y0)
{
  Check_Type(y0, T_ARRAY);

  rb_iv_set(self,  "@dim",     INT2NUM( RARRAY_LEN(y0) ) );
  rb_iv_set(self,  "@y0",      y0);
  rb_iv_set(self,  "@rtol",    rb_float_new(1.0e-9) ); 
  rb_iv_set(self,  "@atol",    rb_float_new(1.0e-9) ); 
  rb_iv_set(self,  "@hmax",    rb_float_new(0.0) ); 
  rb_iv_set(self,  "@inistep", rb_float_new(1.0e-4) ); 
  rb_iv_set(self,  "@monitor", INT2NUM(0) ); 
  rb_iv_set(self,  "@dense",   INT2NUM(0) ); 

  return self;
}

 
VALUE limex_new(VALUE klass, VALUE dim)
{
  VALUE argv[1];

  argv[0] = rb_ary_new2(NUM2LONG(dim));

  rb_obj_call_init(klass, 1, argv);

  return Qnil;
}


void Init_Limex() 
{
  cLimex = rb_define_class("Limex", rb_cObject);
  /* rb_define_singleton_method(cLimex, "new", limex_new, 1); */
  rb_define_method(cLimex, "initialize", limex_init, 1);
  rb_define_method(cLimex, "run", limex_run, 1);
  rb_define_method(cLimex, "srun", limex_srun, 2);
  rb_define_method(cLimex, "y0", limex_y0, 0);
  rb_define_method(cLimex, "y0=", limex_y0_eq, 1);
  rb_define_method(cLimex, "rtol=", limex_rtol_eq, 1);
  rb_define_method(cLimex, "atol=", limex_atol_eq, 1);
  rb_define_method(cLimex, "hmax", limex_hmax, 0);
  rb_define_method(cLimex, "hmax=", limex_hmax_eq, 1);
  rb_define_method(cLimex, "inistep", limex_inistep, 0);
  rb_define_method(cLimex, "inistep=", limex_inistep_eq, 1);
  rb_define_method(cLimex, "monitor", limex_monitor, 0);
  rb_define_method(cLimex, "monitor=", limex_monitor_eq, 1);
  rb_define_method(cLimex, "debug", limex_monitor, 0);
  rb_define_method(cLimex, "debug=", limex_monitor_eq, 1);
  rb_define_method(cLimex, "dense", limex_dense, 0);
  rb_define_method(cLimex, "dense=", limex_dense_eq, 1);
  rb_define_method(cLimex, "tolerance", limex_tolerance, 0);
  rb_define_method(cLimex, "interval", limex_interval, 0);
  rb_define_method(cLimex, "solution", limex_solution, 0);
  rb_define_method(cLimex, "steps", limex_steps, 0);
  rb_define_method(cLimex, "ifail", limex_ifail, 0);
  rb_define_method(cLimex, "status", limex_status, 0);
}

