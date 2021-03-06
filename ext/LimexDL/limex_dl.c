#include <ruby.h>
#include <math.h>
// #include <string.h>  /* for strncmp() */
#include "LIMEX4_3A/LIMEX4_3A.h"
#include "Model_ODE/ydot_LIMEX.h"


VALUE cLimexDL;

/*
static 
void limex_dl_free(void *p)
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

/*
static 
void sfcn(int *n, int *nz,
          double *t, double *y, double *dy,
          double *B, int *ir, int *ic,
          int *info)
{
  *info = -987;
}
*/


static
VALUE limex_dl_srun(VALUE self, VALUE tspan, VALUE pidx)
{
  VALUE arr;

  int     j, k, n, nDAE, nTp, nPdx;
  double  t0, T;
  double  *tp, *z, *dz;
  double  rtol, atol;
  double  h, hMax, hInit;
  int     debug, nDense;
  int     iOpt[32];
  double  rOpt[5];
  int     *iPdx, *iPos, iFail[3];
  int     kOrder;
  double  dense[MAX_NO_EQNS*(2 + MAX_ROW_TAB*(MAX_ROW_TAB+1)/2)];
  double  t1, t2;
  double  *y, *rTol, *aTol;

  VALUE y0;
  VALUE steps;
  VALUE solution;
  VALUE tmparr;

  Check_Type(tspan, T_ARRAY);
  Check_Type(pidx, T_ARRAY);

  // double const EPMACH = 2.0*NUM2DBL( rb_const_get(rb_cFloat, rb_intern("EPSILON")) );

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
  iPdx =    (int *) ALLOCA_N(int, nPdx);

  for (k = 0; k < nPdx; ++k)
  {
    iPdx[k] = NUM2INT( rb_ary_entry(pidx, (long)k) );
  }

  init_var_( &nPdx, iPdx );

  //

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
        slimdherm_( &nDAE, &n, ydot_slimex_, 0, &t0, &T, z, dz, 
                    rTol, aTol, &h, iOpt, rOpt, iPos, 
                    iFail, &kOrder, dense, &t1, &t2
                  );

        // fprintf(stderr,"!!! limex_dl_srun:\n");
        // fprintf(stderr,"!!! k=%d, nDAE=%d, n=%d, t0=%f, t1=%f, t2=%f, tp[k]=%f\n",
        //                 k,nDAE,n,t0,t1,t2,tp[k]);

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
                // if ( !(j < nDAE) && (fabs(z[j]) < EPMACH) )
                // {
                //    z[j] = 0.0;
                // }
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
                // if ( !(j < nDAE) && (fabs(y[j]) < EPMACH) )
                // {
                //    y[j] = 0.0;
                // }
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
      slimdherm_( &nDAE, &n, ydot_slimex_, 0, &t0, &T, z, dz, 
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
          // if ( !(j < nDAE) && (fabs(z[j]) < EPMACH) )
          // {
          //   z[j] = 0.0;
          // }
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

  arr = rb_ary_new();
  rb_ary_push( arr, INT2NUM(iFail[0]) );
  rb_ary_push( arr, INT2NUM(iFail[1]) );
  rb_ary_push( arr, INT2NUM(iFail[2]) );

  rb_iv_set(self, "@ifail", arr);

  return arr;
}


static
VALUE limex_dl_run(VALUE self, VALUE tspan)
{
  VALUE arr;

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

  Check_Type(tspan, T_ARRAY);

  nTp = RARRAY_LEN(tspan);

  if ( nTp < 2 )
  {
    rb_raise(rb_eArgError, "tspan array length wrong (must be 2, at least)");
    return Qnil;
  }
  
  // rb_iv_set(self, "@t0", tStart);
  n  = NUM2INT( rb_iv_get(self, "@dim") );
  // t0 = NUM2DBL(tStart);
  // T  = NUM2DBL(tEnd);

/*
  int     nCom, nPar;
  double  *cbBeg = &sbmlvariables_.start;

  nCom = RARRAY_LEN( rb_iv_get(self, "@compart") );
  nPar = RARRAY_LEN( rb_iv_get(self, "@par") );

  cbBeg += (nCom + n);
  fprintf(stderr,"\n\n***\n%s\n", "INSIDE limex_dl_run(): par =");
  for (j = 0; j < nPar; ++j)
  {
    fprintf(stderr, " %.7f ", *cbBeg++);
  }
  fprintf(stderr,"\n\n");
*/

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
  T  = tp[nTp-1];

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
//        *single step w/hermite* LIMDHERM loop (iOpt[31]==1, nTp > 2)
  //  or  *single step*           LIMEX    loop (iOpt[11]==1, nTp==2)
  //  or  *dense output*          LIMEX    loop (iOpt[11]==2, iOpt[13]==nDense, *and* nTp==2)

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
        limdherm_( &n, ydot_limex_, 0, &t0, &T, z, dz,
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
      limex_( &n, ydot_limex_, 0, &t0, &T, z, dz, 
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

  arr = rb_ary_new();
  rb_ary_push( arr, INT2NUM(iFail[0]) );
  rb_ary_push( arr, INT2NUM(iFail[1]) );
  rb_ary_push( arr, INT2NUM(iFail[2]) );

  rb_iv_set(self, "@ifail", arr);

  return arr;
}


static
VALUE limex_dl_interval(VALUE self)
{
  VALUE arr;

  arr = rb_ary_new();
  rb_ary_push( arr, rb_iv_get(self, "@t0") );
  rb_ary_push( arr, rb_iv_get(self, "@T") );

  return arr;
}


static
VALUE limex_dl_steps(VALUE self)
{
  return rb_iv_get(self, "@steps");
}


static
VALUE limex_dl_solution(VALUE self)
{
  return rb_iv_get(self, "@solution");
}


static
VALUE limex_dl_ifail(VALUE self)
{
  return rb_iv_get(self, "@ifail");
}

static
VALUE limex_dl_status(VALUE self)
{
  VALUE       ifail;
  const char  *msg = "n/a";

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
VALUE limex_dl_modelids(VALUE self)
{
  return rb_iv_get(self, "@version");
}


static
VALUE limex_dl_compart(VALUE self)
{
  return rb_iv_get(self, "@compart");
}


static
VALUE limex_dl_y0(VALUE self)
{
  return rb_iv_get(self, "@y0");
}

static
VALUE limex_dl_y0_eq(VALUE self, VALUE y0)
{
  Check_Type(y0, T_ARRAY);

  if ( NUM2INT( rb_iv_get(self, "@dim") ) != RARRAY_LEN(y0) ) 
  {
    rb_raise(rb_eArgError, "array length wrong");
    return Qnil;
  }

  rb_iv_set(self,  "@dim",  INT2NUM( RARRAY_LEN(y0) ) );
  rb_iv_set(self,  "@y0",   y0);

  return y0;
}

static
VALUE limex_dl_yId(VALUE self)
{
  return rb_iv_get(self, "@yId");
}


static
VALUE limex_dl_par(VALUE self)
{
  // fprintf(stderr, "\n\n***\n%s\n", "INSIDE limex_dl_par()");

  return rb_iv_get(self, "@par");
}

static
VALUE limex_dl_par_eq(VALUE self, VALUE par)
{
  double *p, dummy;
  int    j, k, l, m, nPar;
  VALUE  _par;

  Check_Type(par, T_ARRAY);

  // fprintf(stderr, "\n\n***\n%s\n", "INSIDE limex_dl_par_eq()");

  nPar = RARRAY_LEN(par);
  _par = rb_iv_get(self, "@par");
  
  if ( RARRAY_LEN(_par) != nPar )
  {
    rb_raise(rb_eArgError, "array length wrong"); 
    return Qnil;
  }

  rb_iv_set(self,  "@par",   par);

  p = (double *) ALLOCA_N(double, nPar);
  k = l = m = dummy = 0;
  for (j = 0; j < nPar; ++j)
  {
     p[j] = NUM2DBL( rb_ary_entry(par, (long)j) );
  }
  init_ode_( &dummy,&k,&l, &dummy,&k,&m, p,&k,&nPar ); 

  return par;
}


static
VALUE limex_dl_pId(VALUE self)
{
  return rb_iv_get(self, "@pId");
}


static
VALUE limex_dl_monitor(VALUE self)
{
  return rb_iv_get(self, "@monitor");
}

static
VALUE limex_dl_monitor_eq(VALUE self, VALUE dbg)
{
  rb_iv_set(self,  "@monitor",   dbg);

  return dbg;
}


static
VALUE limex_dl_tolerance(VALUE self)
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
VALUE limex_dl_rtol_eq(VALUE self, VALUE rtol)
{
  Check_Type(rtol, T_FLOAT);

  rb_iv_set(self, "@rtol", rtol);

  return rtol;
}

static
VALUE limex_dl_atol_eq(VALUE self, VALUE atol)
{
  Check_Type(atol, T_FLOAT);

  rb_iv_set(self, "@atol", atol);

  return atol;
}


static
VALUE limex_dl_inistep(VALUE self)
{
  return rb_iv_get(self, "@inistep");
}

static
VALUE limex_dl_inistep_eq(VALUE self, VALUE inistep)
{
  Check_Type(inistep, T_FLOAT);

  rb_iv_set(self, "@inistep", inistep);

  return inistep;
}


static
VALUE limex_dl_hmax(VALUE self)
{
  return rb_iv_get(self, "@hmax");
}

static
VALUE limex_dl_hmax_eq(VALUE self, VALUE hmax)
{
  Check_Type(hmax, T_FLOAT);

  rb_iv_set(self, "@hmax", hmax);

  return hmax;
}


static
VALUE limex_dl_dense(VALUE self)
{
  return rb_iv_get(self, "@dense");
}

static
VALUE limex_dl_dense_eq(VALUE self, VALUE dense)
{
  return rb_iv_set(self, "@dense", dense);
}


static 
VALUE limex_dl_init(VALUE self) /* , VALUE y0) */
{
  /* Check_Type(y0, T_ARRAY); */
  double *cbBeg, dummy;
  char   *idm, *ids, *idp;
  int    j, k, nCom, dim, nPar;
  int    nId = 32;

  VALUE version;
  VALUE com;
  VALUE y0, yId;
  VALUE par, pId;

  dummy = k = nCom = dim = nPar = 0;
 
  /* fprintf(stderr, "\n\n***\n%s\n", "BEFORE init_ode_"); */

  init_ode_( &dummy,&k,&nCom, &dummy,&k,&dim, &dummy,&k,&nPar );

  version = rb_ary_new();
      com = rb_ary_new();
       y0 = rb_ary_new();
      yId = rb_ary_new();
      par = rb_ary_new();
      pId = rb_ary_new();

  cbBeg = &sbmlvariables_.start;

  for (j = 0; j < nCom; ++j)
  {
    rb_ary_push( com, rb_float_new(*cbBeg++) );
  }
  for (j = 0; j < dim; ++j)
  {
    rb_ary_push( y0, rb_float_new(*cbBeg++) );
  }
  for (j = 0; j < nPar; ++j)
  {
    rb_ary_push( par, rb_float_new(*cbBeg++) );
  }

  idm = (char *) ALLOCA_N(char,nId*MAXIDSTRLEN);
  ids = (char *) ALLOCA_N(char,dim*MAXIDSTRLEN);
  idp = (char *) ALLOCA_N(char,nPar*MAXIDSTRLEN);
  /*
  for (j = 0; j < dim; ++j)
  {
    ids[j] = (char *) ALLOCA_N(char, MAXIDSTRLEN);
  }
  for (j = 0; j < nPar; ++j)
  {
    idp[j] = (char *) ALLOCA_N(char, MAXIDSTRLEN);
  }  
  */

  get_model_ids_( idm, &nId, MAXIDSTRLEN );
  get_species_ids_( ids, &dim, MAXIDSTRLEN );
  get_parameter_ids_( idp, &nPar, MAXIDSTRLEN );

  for (j = 0; j < nId; ++j)
  {
    rb_ary_push( version, rb_str_new2(idm+j*MAXIDSTRLEN) );
  }
  for (j = 0; j < dim; ++j)
  {
    rb_ary_push( yId, rb_str_new2(ids+j*MAXIDSTRLEN) );
  }
  for (j = 0; j < nPar; ++j)
  {
    rb_ary_push( pId, rb_str_new2(idp+j*MAXIDSTRLEN) );
  }

  /**/
  /*
  if ( strncmp(idm+4*MAXIDSTRLEN, "with vareq", MAXIDSTRLEN) == 0 )
  {
    set_model_number_directions_();
    k = 0;
    init_var_(&k,NULL);
  }
  */
  /**/

  rb_iv_set(self,  "@dim",     INT2NUM( RARRAY_LEN(y0) ) );
  rb_iv_set(self,  "@version", version);
  rb_iv_set(self,  "@compart", com);
  rb_iv_set(self,  "@y0",      y0);
  rb_iv_set(self,  "@yId",     yId);
  rb_iv_set(self,  "@par",     par);
  rb_iv_set(self,  "@pId",     pId);
  rb_iv_set(self,  "@rtol",    rb_float_new(1.0e-9) ); 
  rb_iv_set(self,  "@atol",    rb_float_new(1.0e-9) ); 
  rb_iv_set(self,  "@hmax",    rb_float_new(0.0) ); 
  rb_iv_set(self,  "@inistep", rb_float_new(1.0e-4) ); 
  rb_iv_set(self,  "@monitor", INT2NUM(0) ); 
  rb_iv_set(self,  "@dense",   INT2NUM(0) ); 

  return self;
}

 
VALUE limex_dl_new(VALUE klass, VALUE dim)
{
  VALUE argv[1];

  argv[0] = rb_ary_new2(NUM2LONG(dim));

  rb_obj_call_init(klass, 1, argv);

  return Qnil;
}


void Init_LimexDL() 
{
  cLimexDL = rb_define_class("LimexDL", rb_cObject);
  /* rb_define_singleton_method(cLimexDL, "new", limex_dl_new, 1); */
  rb_define_method(cLimexDL, "initialize", limex_dl_init, 0);
  rb_define_method(cLimexDL, "version", limex_dl_modelids, 0);
  rb_define_method(cLimexDL, "run", limex_dl_run, 1);
  rb_define_method(cLimexDL, "srun", limex_dl_srun, 2);
  rb_define_method(cLimexDL, "compart", limex_dl_compart, 0);
  rb_define_method(cLimexDL, "y0", limex_dl_y0, 0);
  rb_define_method(cLimexDL, "y0=", limex_dl_y0_eq, 1);
  rb_define_method(cLimexDL, "yId", limex_dl_yId, 0);
  rb_define_method(cLimexDL, "par", limex_dl_par, 0);
  rb_define_method(cLimexDL, "par=", limex_dl_par_eq, 1);
  rb_define_method(cLimexDL, "pId", limex_dl_pId, 0);
  rb_define_method(cLimexDL, "rtol=", limex_dl_rtol_eq, 1);
  rb_define_method(cLimexDL, "atol=", limex_dl_atol_eq, 1);
  rb_define_method(cLimexDL, "hmax", limex_dl_hmax, 0);
  rb_define_method(cLimexDL, "hmax=", limex_dl_hmax_eq, 1);
  rb_define_method(cLimexDL, "inistep", limex_dl_inistep, 0);
  rb_define_method(cLimexDL, "inistep=", limex_dl_inistep_eq, 1);
  rb_define_method(cLimexDL, "monitor", limex_dl_monitor, 0);
  rb_define_method(cLimexDL, "monitor=", limex_dl_monitor_eq, 1);
  rb_define_method(cLimexDL, "debug", limex_dl_monitor, 0);
  rb_define_method(cLimexDL, "debug=", limex_dl_monitor_eq, 1);
  rb_define_method(cLimexDL, "dense", limex_dl_dense, 0);
  rb_define_method(cLimexDL, "dense=", limex_dl_dense_eq, 1);
  rb_define_method(cLimexDL, "tolerance", limex_dl_tolerance, 0);
  rb_define_method(cLimexDL, "interval", limex_dl_interval, 0);
  rb_define_method(cLimexDL, "solution", limex_dl_solution, 0);
  rb_define_method(cLimexDL, "steps", limex_dl_steps, 0);
  rb_define_method(cLimexDL, "ifail", limex_dl_ifail, 0);
  rb_define_method(cLimexDL, "status", limex_dl_status, 0);
}

