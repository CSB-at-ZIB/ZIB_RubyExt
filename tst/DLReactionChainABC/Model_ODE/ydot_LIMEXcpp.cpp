/*
c-----
c SBML Model : ChainABC                                               
c              ~~~~~~~~
c       Date : Thu Jul 30 13:24:00 2015
c              
c     Author : automated transcription by 'sbml2adolc'
c              
c Copyright (C) Zuse Institute Berlin, CSB Group
c-----
*/
/*
   ///
c  com[  0]  :   default                                   
   ///
c
   ///
c  spe[  0]  :   A                                         
c  spe[  1]  :   B                                         
c  spe[  2]  :   C                                         
   ///
c
   ///
c  par[  0]  :   global_k1                                 
c  par[  1]  :   global_k_1                                
c  par[  2]  :   global_k2                                 
c  par[  3]  :   global_k_2                                
   ///
*/
#include <cmath>
#include <cstring>
#include <algorithm>
// do not even think to #include "ydot_LIMEXcpp.h"
/**/
#define ADOLC_TAPELESS
#define NUM_SPE 3
#define NUM_PAR 4
#define NUMBER_DIRECTIONS (NUM_SPE + NUM_PAR)
// ... #include <adolc/adtl.h>
#include <adolc/adouble.h>
typedef adtl::adouble adouble;
ADOLC_TAPELESS_UNIQUE_INTERNALS;
/**/
#define MAXIDSTRLEN 64
/**/
using namespace std;
//=======================================================================
extern "C" {
//=======================================================================
struct
{
      double  com[1];
      double  spe[3];
      double  par[4];
      double  rul[0];
      double  rea[4];
      int     eve[0];
/*
      adouble adcom[1];
      adouble adspe[3];
      adouble adpar[4];
      adouble adrul[0];
      adouble adrea[4];
      int     adeve[0];
      int     apidx[4];
      int     npidx;
*/
} sbmlvariables_;
//=======================================================================
struct {
      adouble com[1];
      adouble spe[3];
      adouble par[4];
      adouble rul[0];
      adouble rea[4];
      int     eve[0];
      int     apidx[4];
      int     npidx;
} advariables_;
//=======================================================================
#define  zero  0.0e0
#define  one   1.0e0
#define  pi    3.141592653589793238462643383276e0
#define  eul   2.718281828459045235360287471352e0
#define  com   PRE(variables_).com
#define  spe   PRE(variables_).spe
#define  par   PRE(variables_).par
#define  rul   PRE(variables_).rul
#define  rea   PRE(variables_).rea
#define  eve   PRE(variables_).eve
#define  apidx advariables_.apidx
#define  npidx advariables_.npidx
//=======================================================================
void compute_fy( int, int, double, double*, double*, double*, int* );
void compute_fp( int, int, double, double*, double*, int* );
void ydotAD ( int, double, adouble*, adouble*, int* );
void check_eventsAD ( double, int, int* ); 
void ydot ( int, double, double*, double*, int* );
void check_events ( double, int, int* ); 
//=======================================================================
};
//=======================================================================
 
//-----------------------------------------------------------------------
 
//=======================================================================
extern "C" {
//=======================================================================
void ydot_slimex_ ( int* n, int* nz, double* t, double* y_state, double* dy,
                    double* b, int* ir, int* ic, int* info )
{
   int    nspe = 3;
   int    npar = npidx;
   /**/
   *nz = *n;
   for (int j = 0; j < *nz; ++j)
   {
      ir[j] =
      ic[j] = j+1;
       b[j] = one;
   }
   /**/
   if ( *n == nspe )
   {
 
      ydot ( *n, *t, y_state, dy, info );
      return;
   }
   /**/
   double fp[npar*nspe];
   double fy[nspe*nspe];
   double  s[*n];
   /**/
   for (int j = 0; j < *n; ++j)
   {
     s[j] = y_state[j];
   }
   /**/
   compute_fp( nspe, npar, *t, y_state, fp, info );
   compute_fy( nspe, npar, *t, y_state, dy, fy, info );
   /**/
   for (int idx = 0; idx < npar; ++idx)
   {
      for (int j = 0; j < nspe; ++j)
      {
          double sum = 0.0;
          if ( apidx[idx] > 0 ) 
          {
              sum = fp[j + nspe*idx]; /* fp[j][ apidx[idx]-1 ]; */
          }
          for (int nu = 0; nu < nspe; ++nu)
          {
              sum += fy[j + nspe*nu] * s[nu + nspe*(idx+1)];
          } 
          dy[j + nspe*(idx+1)] = sum;
      }
   }
}
//=======================================================================
#define PRE(x) ad ## x
void compute_fy( int nspe, int npar, double t, double* y_state, 
                 double* dy, double* fy, int* info )
{
   /// adtl::setNumDir(nspe);
   adouble ay[nspe], ady[nspe];
   /**/
   for (int j = 0; j < nspe; ++j)
   {
      ay[j] = y_state[j];
      ay[j].setADValue(j,1);
   }
   for (int k = 0; k < npar; ++k)
   {
      if ( apidx[k] > 0 )
      {  
         par[ apidx[k]-1 ].setADValue(nspe+k,0);
      }
   }
   /**/
   ydotAD ( nspe, t, ay, ady, info );
   /**/
   for (int k = 0; k < nspe; ++k)
   {
      dy[k] = ady[k].getValue();
      /**/
      for (int j = 0; j < nspe; ++j)
      {
         fy[ j + nspe*k ] = ady[j].getADValue(k);
      }
   }
}
//-----------------------------------------------------------------------
void compute_fp( int nspe, int npar, double t, double* y_state, 
                 double* fp, int* info )
{
   /// adtl::setNumDir(npar);
   adouble ay[nspe], ady[nspe];
   /**/
   for (int j = 0; j < nspe; ++j)
   {
      ay[j] = y_state[j];
   }
   for (int k = 0; k < npar; ++k)
   {
      if ( apidx[k] > 0 )
      {  
         par[ apidx[k]-1 ].setADValue(nspe+k,1);
      }
   }
   /**/
   ydotAD ( nspe, t, ay, ady, info );
   /**/
   for (int k = 0; k < npar; ++k)
   {
      for (int j = 0; j < nspe; ++j)
      {
         fp[ j + nspe*k ] = ady[j].getADValue(nspe+k);
      }
   }
}
#undef PRE
//======================================================================= 
#define PRE(x) ad ## x
void ydotAD ( int n, double t, adouble* y, adouble* dy, int* info)
{
/*
 
 
*/
   *info = 0;
   for (int j = 0; j < n; ++j)
   {
      spe[j] = y[j];
       dy[j] = zero;
   }
//-----------------------------------------------------------------------
   check_eventsAD (t, 0, eve);
//-----------------------------------------------------------------------
 
//-----------------------------------------------------------------------
   check_eventsAD (t, 0, eve);
//-----------------------------------------------------------------------
      rea[0] = par[0] * spe[0];
      rea[1] = par[1] * spe[1];
      rea[2] = par[2] * spe[1];
      rea[3] = par[3] * spe[2];
 
//-----------------------------------------------------------------------
   check_eventsAD (t, 0, eve);
//-----------------------------------------------------------------------
      dy[0] = ( - rea[0] + rea[1] ) / com[0];
      dy[1] = ( + rea[0] - rea[1] - rea[2] + rea[3] ) / com[0];
      dy[2] = ( + rea[2] - rea[3] ) / com[0];
 
}
//======================================================================= 
void check_eventsAD ( double t, int n, int* ev)
{
/*
 
 
*/
 
}
#undef PRE
//======================================================================= 
// #define PRE(x) sbml ## x
#undef com
#undef spe
#undef par
void init_var_ ( int* nidx, int* pidx )
{
   int ncom = 1;
   int nspe = 3;
   int npar = 4;
   /**/
   for (int j = 0; j < ncom; ++j)
   {
      advariables_.com[j] = sbmlvariables_.com[j];
   }
   /**/
   for (int j = 0; j < nspe; ++j)
   {
      advariables_.spe[j] = sbmlvariables_.spe[j];
   }
   /**/
   for (int j = 0; j < npar; ++j)
   {
      advariables_.par[j] = sbmlvariables_.par[j];
   }
   /**/
   int midx = std::min(npar,*nidx);
   for (int idx = 0; idx < midx; ++idx)
   {
      apidx[idx] = pidx[idx];
   }
   npidx = midx;
   /**/
   // adtl::setNumDir(nspe + midx);
}
#define com PRE(variables_).com
#define spe PRE(variables_).spe
#define par PRE(variables_).par
// #undef PRE
//======================================================================= 
void ydot_limex_ ( int* n, int* nz, double* t, double* y_state, double* dy, 
                   double* b, int* ir, int* ic, int* info )
{
   *nz = *n;
   for (int j = 0; j < *nz; ++j)
   {
      ir[j] =
      ic[j] = j+1;
       b[j] = one;
   }
 
   ydot ( *n, *t, y_state, dy, info );
}
//=======================================================================
#define PRE(x) sbml ## x
void ydot ( int n, double t, double* y, double* dy, int* info )
{
/*
 
 
*/
   *info = 0;
   for (int j = 0; j < n; ++j)
   {
      spe[j] = y[j];
       dy[j] = zero;
   }
//-----------------------------------------------------------------------
   check_events (t, 0, eve);
//-----------------------------------------------------------------------
 
//-----------------------------------------------------------------------
   check_events (t, 0, eve);
//-----------------------------------------------------------------------
      rea[0] = par[0] * spe[0];
      rea[1] = par[1] * spe[1];
      rea[2] = par[2] * spe[1];
      rea[3] = par[3] * spe[2];
 
//-----------------------------------------------------------------------
   check_events (t, 0, eve);
//-----------------------------------------------------------------------
      dy[0] = ( - rea[0] + rea[1] ) / com[0];
      dy[1] = ( + rea[0] - rea[1] - rea[2] + rea[3] ) / com[0];
      dy[2] = ( + rea[2] - rea[3] ) / com[0];
 
}
//=======================================================================
void check_events ( double t, int n, int* ev )
{
/*
 
 
*/
 
}
#undef PRE
//=======================================================================
#define PRE(x) sbml ## x
void init_ode_ ( double* c, int* ic, int* lc ,  
                 double* s, int* is, int* ls ,  
                 double* p, int* ip, int* lp )
{
   com[  0] = 1.000000e+00;
 
   spe[  0] = 1.000000e+00;
   spe[  1] = 0.000000e+00;
   spe[  2] = 0.000000e+00;
 
   par[  0] = 2.000000e+00;
   par[  1] = 3.000000e-03;
   par[  2] = 1.000000e+00;
   par[  3] = 2.000000e-03;
 
 
//-----------------------------------------------------------------------
   int mc = std::min(1,*lc);
   if ( ic[0] == 0 )
   {
      for (int j = 0; j < mc; ++j)
      {
         com[ j ] = c[j];
      }
   }
   else
   {
      for (int j = 0; j < mc; ++j)
      {
         com[ ic[j]-1 ] = c[j];
      }
   }
   if ( *lc == 0 )
   {
      *lc = 1;
   }
//-----------------------------------------------------------------------
   int ms = std::min(3,*ls);
   if ( is[0] == 0 )
   {
      for (int j = 0; j < ms; ++j)
      {
            spe[ j ] = s[j];
      }
   }
   else
   {
      for (int j = 0; j <  ms; ++j)
      {
         spe[ is[j]-1 ] = s[j];
      }
   }
   if ( *ls == 0 )
   {
      *ls = 3;
   }
 
//-----------------------------------------------------------------------
   int mp = std::min(4,*lp);
   if ( ip[0] == 0 )
   {
      for (int j = 0; j < mp; ++j)
      {
         par[ j ] = p[j];
      }
   }
   else
   {
      for (int j = 0; j < mp; ++j)
      {
         par[ ip[j]-1 ] = p[j];
      }
   }
   if ( *lp == 0 )
   {
      *lp = 4;
   }
}
#undef PRE
//=======================================================================
void get_compartment_ids_ ( char (*idc)[MAXIDSTRLEN], int* ncom, int _len)
{
   *ncom = 1;
   ///
   strncpy(idc[  0],"default\0",_len);
   ///
}
//=======================================================================
void get_species_ids_ ( char (*ids)[MAXIDSTRLEN], int* nspe, int _len)
{
   *nspe = 3;
   for (int j = 0; j < *nspe; ++j)
   {
      strncpy(ids[j],"\0",_len);
   }
   ///
   strncpy(ids[  0],"A\0",_len);
   strncpy(ids[  1],"B\0",_len);
   strncpy(ids[  2],"C\0",_len);
   ///
}
//=======================================================================
void get_parameter_ids_ ( char (*idp)[MAXIDSTRLEN], int* npar, int _len)
{
   *npar = 4;
   for (int j = 0; j < *npar; ++j)
   {
      strncpy(idp[j],"\0",_len);
   }
   ///
   strncpy(idp[  0],"global_k1\0",_len);
   strncpy(idp[  1],"global_k_1\0",_len);
   strncpy(idp[  2],"global_k2\0",_len);
   strncpy(idp[  3],"global_k_2\0",_len);
   ///
}
//=======================================================================
void get_model_ids_ ( char (*idm)[MAXIDSTRLEN], int* nid, int _len)
{
   *nid = 5;
   for (int j = 0; j < *nid; ++j)
   {
      strncpy(idm[j],"\0",_len);
   }
   ///
   strncpy(idm[  0],"ChainABC\0",_len);
   strncpy(idm[  1],"Thu Jul 30 13:24:00 2015\0",_len);
   strncpy(idm[  2],"001438255440\0",_len);
   strncpy(idm[  3],"ChainABC.xml\0",_len);
   strncpy(idm[  4],"with vareq\0",_len);
   ///
}
//=======================================================================
};
//=======================================================================
 
 
//=======================================================================

