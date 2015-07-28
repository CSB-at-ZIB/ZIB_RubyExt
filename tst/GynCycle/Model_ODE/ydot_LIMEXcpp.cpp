/*
c-----
c SBML Model : PAEON_V2                                               
c              ~~~~~~~~
c       Date : Tue Jul 28 19:03:47 2015
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
c  spe[  0]  :   LH_pit                                    
c  spe[  1]  :   LH_blood                                  
c  spe[  2]  :   R_LH                                      
c  spe[  3]  :   LH_R                                      
c  spe[  4]  :   R_LH_des                                  
c  spe[  5]  :   FSH_pit                                   
c  spe[  6]  :   FSH_blood                                 
c  spe[  7]  :   R_FSH                                     
c  spe[  8]  :   FSH_R                                     
c  spe[  9]  :   R_FSH_des                                 
c  spe[ 10]  :   S_foll                                    
c  spe[ 11]  :   AF1                                       
c  spe[ 12]  :   AF2                                       
c  spe[ 13]  :   AF3                                       
c  spe[ 14]  :   AF4                                       
c  spe[ 15]  :   PrF                                       
c  spe[ 16]  :   OvF                                       
c  spe[ 17]  :   Sc1                                       
c  spe[ 18]  :   Sc2                                       
c  spe[ 19]  :   Lut1                                      
c  spe[ 20]  :   Lut2                                      
c  spe[ 21]  :   Lut3                                      
c  spe[ 22]  :   Lut4                                      
c  spe[ 23]  :   E2                                        
c  spe[ 24]  :   P4                                        
c  spe[ 25]  :   InhA                                      
c  spe[ 26]  :   InhB                                      
c  spe[ 27]  :   InhA_tau                                  
c  spe[ 28]  :   GnRH_G                                    
c  spe[ 29]  :   R_GnRH_a                                  
c  spe[ 30]  :   R_GnRH_i                                  
c  spe[ 31]  :   GnRH_R_a                                  
c  spe[ 32]  :   GnRH_R_i                                  
   ///
c
   ///
c  par[  0]  :   global_p_001_001                          
c  par[  1]  :   global_p_001_002                          
c  par[  2]  :   global_p_001_003                          
c  par[  3]  :   global_p_001_004                          
c  par[  4]  :   global_p_001_005                          
c  par[  5]  :   global_p_001_006                          
c  par[  6]  :   global_p_001_007                          
c  par[  7]  :   global_p_001_008                          
c  par[  8]  :   global_p_001_009                          
c  par[  9]  :   global_p_001_010                          
c  par[ 10]  :   global_p_002_001                          
c  par[ 11]  :   global_p_002_002                          
c  par[ 12]  :   global_p_002_003                          
c  par[ 13]  :   global_p_003_001                          
c  par[ 14]  :   global_p_004_001                          
c  par[ 15]  :   global_p_006_001                          
c  par[ 16]  :   global_p_006_002                          
c  par[ 17]  :   global_p_006_003                          
c  par[ 18]  :   global_p_006_004                          
c  par[ 19]  :   global_p_006_005                          
c  par[ 20]  :   global_p_006_006                          
c  par[ 21]  :   global_p_006_007                          
c  par[ 22]  :   global_p_006_008                          
c  par[ 23]  :   global_p_006_009                          
c  par[ 24]  :   global_p_006_010                          
c  par[ 25]  :   global_p_006_011                          
c  par[ 26]  :   global_p_007_001                          
c  par[ 27]  :   global_p_007_002                          
c  par[ 28]  :   global_p_008_001                          
c  par[ 29]  :   global_p_009_001                          
c  par[ 30]  :   global_p_011_001                          
c  par[ 31]  :   global_p_011_002                          
c  par[ 32]  :   global_p_011_003                          
c  par[ 33]  :   global_p_011_004                          
c  par[ 34]  :   global_p_011_005                          
c  par[ 35]  :   global_p_011_006                          
c  par[ 36]  :   global_p_011_007                          
c  par[ 37]  :   global_p_012_001                          
c  par[ 38]  :   global_p_012_002                          
c  par[ 39]  :   global_p_012_003                          
c  par[ 40]  :   global_p_012_004                          
c  par[ 41]  :   global_p_013_001                          
c  par[ 42]  :   global_p_013_002                          
c  par[ 43]  :   global_p_013_003                          
c  par[ 44]  :   global_p_014_001                          
c  par[ 45]  :   global_p_014_002                          
c  par[ 46]  :   global_p_014_003                          
c  par[ 47]  :   global_p_014_004                          
c  par[ 48]  :   global_p_015_001                          
c  par[ 49]  :   global_p_015_002                          
c  par[ 50]  :   global_p_015_003                          
c  par[ 51]  :   global_p_016_001                          
c  par[ 52]  :   global_p_016_002                          
c  par[ 53]  :   global_p_017_001                          
c  par[ 54]  :   global_p_017_002                          
c  par[ 55]  :   global_p_017_003                          
c  par[ 56]  :   global_p_017_004                          
c  par[ 57]  :   global_p_018_001                          
c  par[ 58]  :   global_p_018_002                          
c  par[ 59]  :   global_p_018_003                          
c  par[ 60]  :   global_p_018_004                          
c  par[ 61]  :   global_p_019_001                          
c  par[ 62]  :   global_p_020_001                          
c  par[ 63]  :   global_p_020_002                          
c  par[ 64]  :   global_p_020_003                          
c  par[ 65]  :   global_p_020_004                          
c  par[ 66]  :   global_p_021_001                          
c  par[ 67]  :   global_p_022_001                          
c  par[ 68]  :   global_p_023_001                          
c  par[ 69]  :   global_p_024_001                          
c  par[ 70]  :   global_p_024_002                          
c  par[ 71]  :   global_p_024_003                          
c  par[ 72]  :   global_p_024_004                          
c  par[ 73]  :   global_p_024_005                          
c  par[ 74]  :   global_p_024_006                          
c  par[ 75]  :   global_p_024_007                          
c  par[ 76]  :   global_p_024_008                          
c  par[ 77]  :   global_p_025_001                          
c  par[ 78]  :   global_p_025_002                          
c  par[ 79]  :   global_p_025_003                          
c  par[ 80]  :   global_p_026_001                          
c  par[ 81]  :   global_p_026_002                          
c  par[ 82]  :   global_p_026_003                          
c  par[ 83]  :   global_p_026_004                          
c  par[ 84]  :   global_p_026_005                          
c  par[ 85]  :   global_p_026_006                          
c  par[ 86]  :   global_p_026_007                          
c  par[ 87]  :   global_p_026_008                          
c  par[ 88]  :   global_p_027_001                          
c  par[ 89]  :   global_p_027_002                          
c  par[ 90]  :   global_p_027_003                          
c  par[ 91]  :   global_p_027_004                          
c  par[ 92]  :   global_p_028_001                          
c  par[ 93]  :   global_p_029_001                          
c  par[ 94]  :   global_p_029_002                          
c  par[ 95]  :   global_p_029_003                          
c  par[ 96]  :   global_p_029_004                          
c  par[ 97]  :   global_p_029_005                          
c  par[ 98]  :   global_p_029_006                          
c  par[ 99]  :   global_p_030_001                          
c  par[100]  :   global_p_030_002                          
c  par[101]  :   global_p_030_003                          
c  par[102]  :   global_p_030_004                          
c  par[103]  :   global_p_030_005                          
c  par[104]  :   global_p_031_001                          
c  par[105]  :   global_p_031_002                          
c  par[106]  :   global_p_031_003                          
c  par[107]  :   global_p_032_001                          
c  par[108]  :   global_p_032_002                          
c  par[109]  :   global_p_033_001                          
c  par[110]  :   global_p_033_002                          
c  par[111]  :   global_p_033_003                          
c  par[112]  :   global_p_034_001                          
c  par[113]  :   global_p_034_002                          
c  par[114]  :   global_p_035_001                          
   ///
*/
#include <cmath>
#include <cstring>
#include <algorithm>
// do not even think to #include "ydot_LIMEXcpp.h"
/**/
#define ADOLC_TAPELESS
#define NUM_SPE 33
#define NUM_PAR 115
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
      double  spe[33];
      double  par[115];
      double  rul[2];
      double  rea[57];
      int     eve[0];
      /**/
      adouble adcom[1];
      adouble adspe[33];
      adouble adpar[115];
      adouble adrul[2];
      adouble adrea[57];
      int     adeve[0];
      int     apidx[115];
      int     npidx;
} sbmlvariables_;
//=======================================================================
/*
struct {
      adouble adcom[1];
      adouble adspe[33];
      adouble adpar[115];
      adouble adrul[2];
      adouble adrea[57];
      int     adeve[0];
} advariables_;
*/
//=======================================================================
#define  zero  0.0e0
#define  one   1.0e0
#define  pi    3.141592653589793238462643383276e0
#define  eul   2.718281828459045235360287471352e0
#define  com   sbmlvariables_.PRE(com)
#define  spe   sbmlvariables_.PRE(spe)
#define  par   sbmlvariables_.PRE(par)
#define  rul   sbmlvariables_.PRE(rul)
#define  rea   sbmlvariables_.PRE(rea)
#define  eve   sbmlvariables_.PRE(eve)
#define  apidx sbmlvariables_.apidx
#define  npidx sbmlvariables_.npidx
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
      template <typename T> T fun0 ( T& x1, T& x2, T& x3 );
 
      template <typename T> T fun1 ( T& x1, T& x2, T& x3 );
 
 
//-----------------------------------------------------------------------
 
//=======================================================================
extern "C" {
//=======================================================================
void ydot_slimex_ ( int* n, int* nz, double* t, double* y_state, double* dy,
                    double* b, int* ir, int* ic, int* info )
{
   int    nspe = 33;
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
   for (int idx = 0; idx < npidx; ++idx)
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
      par[k].setADValue(k,0);
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
   adouble fun0(...);
   adouble fun1(...);
 
 
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
      rul[0] = par[93] * fun1(spe[24], par[94], par[95]) *
                (1.000000e+00 + par[96] * fun0(spe[23], par[97],
                par[98]));
      rul[1] = par[99] * (fun0(spe[23], par[100], par[101]) +
                fun1(spe[23], par[102], par[103]));
 
//-----------------------------------------------------------------------
   check_eventsAD (t, 0, eve);
//-----------------------------------------------------------------------
      rea[0] = (par[0] + par[1] * fun0(spe[23], par[2], par[3])) *
                fun1(spe[24], par[4], par[5]);
      rea[1] = (par[6] + par[7] * fun0(spe[31], par[8], par[9])) *
                spe[0];
      rea[2] = (par[6] + par[7] * fun0(spe[31], par[8], par[9])) *
                spe[0] / par[10];
      rea[3] = par[11] * spe[1] * spe[2];
      rea[4] = par[12] * spe[1];
      rea[5] = par[13] * spe[4];
      rea[6] = par[14] * spe[3];
      rea[7] = par[15] / (1.000000e+00 + pow(spe[27] / par[16],
                par[17]) + pow(spe[26] / par[18], par[19])) *
                fun1(rul[0], par[20], par[21]);
      rea[8] = (par[22] + par[23] * fun0(spe[31], par[24], par[25]))
                * spe[5];
      rea[9] = (par[22] + par[23] * fun0(spe[31], par[24], par[25]))
                * spe[5] / par[10];
      rea[10] = par[26] * spe[7] * spe[6];
      rea[11] = par[27] * spe[6];
      rea[12] = par[28] * spe[9];
      rea[13] = par[29] * spe[8];
      rea[14] = par[30] * fun0(spe[6], par[31], par[32]);
      rea[15] = par[33] * fun0(spe[24], par[34], par[35]) * spe[10];
      rea[16] = par[37] * fun0(spe[8], par[38], par[39]);
      rea[17] = par[40] * spe[8] * spe[11];
      rea[18] = par[41] * pow(spe[3] / par[42], par[43]) * spe[10] *
                spe[12];
      rea[19] = par[44] * spe[8] * spe[13] * (1.000000e+00 - spe[13]
                / par[45]);
      rea[20] = par[46] * pow(spe[3] / par[42], par[47]) * spe[10] *
                spe[13];
      rea[21] = par[48] * pow(spe[3] / par[42], par[49]) * spe[14] *
                (1.000000e+00 - spe[14] / par[45]);
      rea[22] = par[50] * (spe[3] / par[42]) * spe[10] * spe[14];
      rea[23] = par[51] * pow(spe[3] / par[42], par[52]) * spe[10] *
                spe[15];
      rea[24] = par[53] * pow(spe[3] / par[42], par[52]) * spe[10] *
                fun0(spe[15], par[54], par[55]);
      rea[25] = par[56] * spe[16];
      rea[26] = par[57] * fun0(spe[16], par[58], par[59]);
      rea[27] = par[60] * spe[17];
      rea[28] = par[61] * spe[18];
      rea[29] = par[62] * spe[19];
      rea[30] = par[62] * par[63] * fun0(spe[31], par[64], par[65])
                * spe[19];
      rea[31] = par[66] * spe[20];
      rea[32] = par[66] * par[63] * fun0(spe[31], par[64], par[65])
                * spe[20];
      rea[33] = par[67] * spe[21];
      rea[34] = par[67] * par[63] * fun0(spe[31], par[64], par[65])
                * spe[21];
      rea[35] = par[68] * (1.000000e+00 + par[63] * fun0(spe[31],
                par[64], par[65])) * spe[22];
      rea[36] = par[69] + par[70] * spe[12] + par[71] * spe[1] *
                spe[13] + par[72] * spe[14] + par[73] * spe[1] *
                spe[15] + par[74] * spe[19] + par[75] * spe[22];
      rea[37] = par[76] * spe[23];
      rea[38] = par[77] + par[78] * spe[22];
      rea[39] = par[79] * spe[24];
      rea[40] = par[80] + par[81] * spe[15] + par[82] * spe[17] +
                par[83] * spe[19] + par[84] * spe[20] + par[85] *
                spe[21] + par[86] * spe[22];
      rea[41] = par[87] * spe[25];
      rea[42] = par[88] + par[89] * spe[12] + par[90] * spe[18];
      rea[43] = par[91] * spe[26];
      rea[44] = par[92] * spe[27];
      rea[45] = rul[1] * rul[0];
      rea[46] = par[104] * spe[29] * spe[28];
      rea[47] = par[105] * spe[31];
      rea[48] = par[106] * spe[28];
      rea[49] = par[107] * spe[29];
      rea[50] = par[108] * spe[30];
      rea[51] = par[109];
      rea[52] = par[110] * spe[32];
      rea[53] = par[111] * spe[30];
      rea[54] = par[112] * spe[31];
      rea[55] = par[113] * spe[32];
      rea[56] = par[114] * spe[32];
 
//-----------------------------------------------------------------------
   check_eventsAD (t, 0, eve);
//-----------------------------------------------------------------------
      dy[0] = ( + rea[0] - rea[1] ) / com[0];
      dy[1] = ( + rea[2] - rea[3] - rea[4] ) / com[0];
      dy[2] = ( - rea[3] + rea[5] ) / com[0];
      dy[3] = ( + rea[3] - rea[6] ) / com[0];
      dy[4] = ( - rea[5] + rea[6] ) / com[0];
      dy[5] = ( + rea[7] - rea[8] ) / com[0];
      dy[6] = ( + rea[9] - rea[10] - rea[11] ) / com[0];
      dy[7] = ( - rea[10] + rea[12] ) / com[0];
      dy[8] = ( + rea[10] - rea[13] ) / com[0];
      dy[9] = ( - rea[12] + rea[13] ) / com[0];
      dy[10] = ( + rea[14] - rea[15] ) / com[0];
      dy[11] = ( + rea[16] - rea[17] ) / com[0];
      dy[12] = ( + rea[17] - rea[18] ) / com[0];
      dy[13] = ( + rea[18] + rea[19] - rea[20] ) / com[0];
      dy[14] = ( + rea[20] + rea[21] - rea[22] ) / com[0];
      dy[15] = ( + rea[22] - rea[23] ) / com[0];
      dy[16] = ( + rea[24] - rea[25] ) / com[0];
      dy[17] = ( + rea[26] - rea[27] ) / com[0];
      dy[18] = ( + rea[27] - rea[28] ) / com[0];
      dy[19] = ( + rea[28] - rea[29] - rea[30] ) / com[0];
      dy[20] = ( + rea[29] - rea[31] - rea[32] ) / com[0];
      dy[21] = ( + rea[31] - rea[33] - rea[34] ) / com[0];
      dy[22] = ( + rea[33] - rea[35] ) / com[0];
      dy[23] = ( + rea[36] - rea[37] ) / com[0];
      dy[24] = ( + rea[38] - rea[39] ) / com[0];
      dy[25] = ( + rea[40] - rea[41] ) / com[0];
      dy[26] = ( + rea[42] - rea[43] ) / com[0];
      dy[27] = ( + rea[41] - rea[44] ) / com[0];
      dy[28] = ( + rea[45] - rea[46] + rea[47] - rea[48] ) / com[0];
      dy[29] = ( - rea[46] + rea[47] - rea[49] + rea[50] ) / com[0];
      dy[30] = ( + rea[49] - rea[50] + rea[51] + rea[52] - rea[53] )
                / com[0];
      dy[31] = ( + rea[46] - rea[47] - rea[54] + rea[55] ) / com[0];
      dy[32] = ( - rea[52] + rea[54] - rea[55] - rea[56] ) / com[0];
 
}
//======================================================================= 
void check_eventsAD ( double t, int n, int* ev)
{
/*
   adouble fun0(...);
   adouble fun1(...);
 
 
*/
 
}
#undef PRE
//======================================================================= 
#define PRE(x) x
void init_var_ ( int* nidx, int* pidx )
{
   int ncom = 1;
   int nspe = 33;
   int npar = 115;
   /**/
   for (int j = 0; j < ncom; ++j)
   {
      sbmlvariables_.adcom[j] = com[j];
   }
   /**/
   for (int j = 0; j < nspe; ++j)
   {
      sbmlvariables_.adspe[j] = spe[j];
   }
   /**/
   for (int j = 0; j < npar; ++j)
   {
      sbmlvariables_.adpar[j] = par[j];
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
#undef PRE
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
#define PRE(x) x
void ydot ( int n, double t, double* y, double* dy, int* info )
{
/*
   T fun0(...);
   T fun1(...);
 
 
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
      rul[0] = par[93] * fun1(spe[24], par[94], par[95]) *
                (1.000000e+00 + par[96] * fun0(spe[23], par[97],
                par[98]));
      rul[1] = par[99] * (fun0(spe[23], par[100], par[101]) +
                fun1(spe[23], par[102], par[103]));
 
//-----------------------------------------------------------------------
   check_events (t, 0, eve);
//-----------------------------------------------------------------------
      rea[0] = (par[0] + par[1] * fun0(spe[23], par[2], par[3])) *
                fun1(spe[24], par[4], par[5]);
      rea[1] = (par[6] + par[7] * fun0(spe[31], par[8], par[9])) *
                spe[0];
      rea[2] = (par[6] + par[7] * fun0(spe[31], par[8], par[9])) *
                spe[0] / par[10];
      rea[3] = par[11] * spe[1] * spe[2];
      rea[4] = par[12] * spe[1];
      rea[5] = par[13] * spe[4];
      rea[6] = par[14] * spe[3];
      rea[7] = par[15] / (1.000000e+00 + pow(spe[27] / par[16],
                par[17]) + pow(spe[26] / par[18], par[19])) *
                fun1(rul[0], par[20], par[21]);
      rea[8] = (par[22] + par[23] * fun0(spe[31], par[24], par[25]))
                * spe[5];
      rea[9] = (par[22] + par[23] * fun0(spe[31], par[24], par[25]))
                * spe[5] / par[10];
      rea[10] = par[26] * spe[7] * spe[6];
      rea[11] = par[27] * spe[6];
      rea[12] = par[28] * spe[9];
      rea[13] = par[29] * spe[8];
      rea[14] = par[30] * fun0(spe[6], par[31], par[32]);
      rea[15] = par[33] * fun0(spe[24], par[34], par[35]) * spe[10];
      rea[16] = par[37] * fun0(spe[8], par[38], par[39]);
      rea[17] = par[40] * spe[8] * spe[11];
      rea[18] = par[41] * pow(spe[3] / par[42], par[43]) * spe[10] *
                spe[12];
      rea[19] = par[44] * spe[8] * spe[13] * (1.000000e+00 - spe[13]
                / par[45]);
      rea[20] = par[46] * pow(spe[3] / par[42], par[47]) * spe[10] *
                spe[13];
      rea[21] = par[48] * pow(spe[3] / par[42], par[49]) * spe[14] *
                (1.000000e+00 - spe[14] / par[45]);
      rea[22] = par[50] * (spe[3] / par[42]) * spe[10] * spe[14];
      rea[23] = par[51] * pow(spe[3] / par[42], par[52]) * spe[10] *
                spe[15];
      rea[24] = par[53] * pow(spe[3] / par[42], par[52]) * spe[10] *
                fun0(spe[15], par[54], par[55]);
      rea[25] = par[56] * spe[16];
      rea[26] = par[57] * fun0(spe[16], par[58], par[59]);
      rea[27] = par[60] * spe[17];
      rea[28] = par[61] * spe[18];
      rea[29] = par[62] * spe[19];
      rea[30] = par[62] * par[63] * fun0(spe[31], par[64], par[65])
                * spe[19];
      rea[31] = par[66] * spe[20];
      rea[32] = par[66] * par[63] * fun0(spe[31], par[64], par[65])
                * spe[20];
      rea[33] = par[67] * spe[21];
      rea[34] = par[67] * par[63] * fun0(spe[31], par[64], par[65])
                * spe[21];
      rea[35] = par[68] * (1.000000e+00 + par[63] * fun0(spe[31],
                par[64], par[65])) * spe[22];
      rea[36] = par[69] + par[70] * spe[12] + par[71] * spe[1] *
                spe[13] + par[72] * spe[14] + par[73] * spe[1] *
                spe[15] + par[74] * spe[19] + par[75] * spe[22];
      rea[37] = par[76] * spe[23];
      rea[38] = par[77] + par[78] * spe[22];
      rea[39] = par[79] * spe[24];
      rea[40] = par[80] + par[81] * spe[15] + par[82] * spe[17] +
                par[83] * spe[19] + par[84] * spe[20] + par[85] *
                spe[21] + par[86] * spe[22];
      rea[41] = par[87] * spe[25];
      rea[42] = par[88] + par[89] * spe[12] + par[90] * spe[18];
      rea[43] = par[91] * spe[26];
      rea[44] = par[92] * spe[27];
      rea[45] = rul[1] * rul[0];
      rea[46] = par[104] * spe[29] * spe[28];
      rea[47] = par[105] * spe[31];
      rea[48] = par[106] * spe[28];
      rea[49] = par[107] * spe[29];
      rea[50] = par[108] * spe[30];
      rea[51] = par[109];
      rea[52] = par[110] * spe[32];
      rea[53] = par[111] * spe[30];
      rea[54] = par[112] * spe[31];
      rea[55] = par[113] * spe[32];
      rea[56] = par[114] * spe[32];
 
//-----------------------------------------------------------------------
   check_events (t, 0, eve);
//-----------------------------------------------------------------------
      dy[0] = ( + rea[0] - rea[1] ) / com[0];
      dy[1] = ( + rea[2] - rea[3] - rea[4] ) / com[0];
      dy[2] = ( - rea[3] + rea[5] ) / com[0];
      dy[3] = ( + rea[3] - rea[6] ) / com[0];
      dy[4] = ( - rea[5] + rea[6] ) / com[0];
      dy[5] = ( + rea[7] - rea[8] ) / com[0];
      dy[6] = ( + rea[9] - rea[10] - rea[11] ) / com[0];
      dy[7] = ( - rea[10] + rea[12] ) / com[0];
      dy[8] = ( + rea[10] - rea[13] ) / com[0];
      dy[9] = ( - rea[12] + rea[13] ) / com[0];
      dy[10] = ( + rea[14] - rea[15] ) / com[0];
      dy[11] = ( + rea[16] - rea[17] ) / com[0];
      dy[12] = ( + rea[17] - rea[18] ) / com[0];
      dy[13] = ( + rea[18] + rea[19] - rea[20] ) / com[0];
      dy[14] = ( + rea[20] + rea[21] - rea[22] ) / com[0];
      dy[15] = ( + rea[22] - rea[23] ) / com[0];
      dy[16] = ( + rea[24] - rea[25] ) / com[0];
      dy[17] = ( + rea[26] - rea[27] ) / com[0];
      dy[18] = ( + rea[27] - rea[28] ) / com[0];
      dy[19] = ( + rea[28] - rea[29] - rea[30] ) / com[0];
      dy[20] = ( + rea[29] - rea[31] - rea[32] ) / com[0];
      dy[21] = ( + rea[31] - rea[33] - rea[34] ) / com[0];
      dy[22] = ( + rea[33] - rea[35] ) / com[0];
      dy[23] = ( + rea[36] - rea[37] ) / com[0];
      dy[24] = ( + rea[38] - rea[39] ) / com[0];
      dy[25] = ( + rea[40] - rea[41] ) / com[0];
      dy[26] = ( + rea[42] - rea[43] ) / com[0];
      dy[27] = ( + rea[41] - rea[44] ) / com[0];
      dy[28] = ( + rea[45] - rea[46] + rea[47] - rea[48] ) / com[0];
      dy[29] = ( - rea[46] + rea[47] - rea[49] + rea[50] ) / com[0];
      dy[30] = ( + rea[49] - rea[50] + rea[51] + rea[52] - rea[53] )
                / com[0];
      dy[31] = ( + rea[46] - rea[47] - rea[54] + rea[55] ) / com[0];
      dy[32] = ( - rea[52] + rea[54] - rea[55] - rea[56] ) / com[0];
 
}
//=======================================================================
void check_events ( double t, int n, int* ev )
{
/*
   T fun0(...);
   T fun1(...);
 
 
*/
 
}
#undef PRE
//=======================================================================
#define PRE(x) x
void init_ode_ ( double* c, int* ic, int* lc ,  
                 double* s, int* is, int* ls ,  
                 double* p, int* ip, int* lp )
{
   com[  0] = 1.000000e+00;
 
   spe[  0] = 2.596960e+05;
   spe[  1] = 2.710000e+00;
   spe[  2] = 8.398000e+00;
   spe[  3] = 2.660000e-01;
   spe[  4] = 7.100000e-01;
   spe[  5] = 5.973400e+04;
   spe[  6] = 5.280000e+00;
   spe[  7] = 5.903000e+00;
   spe[  8] = 7.960000e-01;
   spe[  9] = 1.800000e+00;
   spe[ 10] = 1.560000e-01;
   spe[ 11] = 2.593000e+00;
   spe[ 12] = 2.314000e+01;
   spe[ 13] = 4.920000e-01;
   spe[ 14] = 1.610000e-05;
   spe[ 15] = 2.449000e-01;
   spe[ 16] = 0.000000e+00;
   spe[ 17] = 1.000000e-08;
   spe[ 18] = 2.000000e-06;
   spe[ 19] = 2.000000e-05;
   spe[ 20] = 3.000000e-04;
   spe[ 21] = 3.000000e-03;
   spe[ 22] = 1.100000e-02;
   spe[ 23] = 4.136000e+01;
   spe[ 24] = 1.978000e+00;
   spe[ 25] = 9.300000e-01;
   spe[ 26] = 6.053000e+01;
   spe[ 27] = 7.114000e+01;
   spe[ 28] = 1.520000e-02;
   spe[ 29] = 9.000000e-03;
   spe[ 30] = 9.600000e-04;
   spe[ 31] = 6.500000e-05;
   spe[ 32] = 5.900000e-05;
 
   par[  0] = 7.309920e+03;
   par[  1] = 7.309920e+03;
   par[  2] = 1.922000e+02;
   par[  3] = 1.000000e+01;
   par[  4] = 2.371000e+00;
   par[  5] = 1.000000e+00;
   par[  6] = 4.760000e-03;
   par[  7] = 1.904000e-01;
   par[  8] = 3.000000e-04;
   par[  9] = 5.000000e+00;
   par[ 10] = 5.000000e+00;
   par[ 11] = 2.143000e+00;
   par[ 12] = 7.485100e+01;
   par[ 13] = 6.894900e+01;
   par[ 14] = 1.833600e+02;
   par[ 15] = 2.213000e+04;
   par[ 16] = 9.581000e+01;
   par[ 17] = 5.000000e+00;
   par[ 18] = 7.000000e+01;
   par[ 19] = 3.000000e+00;
   par[ 20] = 1.000000e+01;
   par[ 21] = 3.000000e+00;
   par[ 22] = 5.700000e-02;
   par[ 23] = 2.720000e-01;
   par[ 24] = 3.000000e-04;
   par[ 25] = 3.000000e+00;
   par[ 26] = 3.529000e+00;
   par[ 27] = 1.142500e+02;
   par[ 28] = 6.102900e+01;
   par[ 29] = 1.383000e+02;
   par[ 30] = 2.190000e-01;
   par[ 31] = 3.000000e+00;
   par[ 32] = 5.000000e+00;
   par[ 33] = 1.343000e+00;
   par[ 34] = 1.235000e+00;
   par[ 35] = 5.000000e+00;
   par[ 36] = 2.000000e-01;
   par[ 37] = 3.663000e+00;
   par[ 38] = 6.080000e-01;
   par[ 39] = 3.000000e+00;
   par[ 40] = 1.221000e+00;
   par[ 41] = 4.882000e+00;
   par[ 42] = 2.726000e+00;
   par[ 43] = 3.689000e+00;
   par[ 44] = 1.220000e-01;
   par[ 45] = 1.000000e+01;
   par[ 46] = 1.220600e+02;
   par[ 47] = 5.000000e+00;
   par[ 48] = 1.220600e+01;
   par[ 49] = 2.000000e+00;
   par[ 50] = 3.327500e+02;
   par[ 51] = 1.220600e+02;
   par[ 52] = 6.000000e+00;
   par[ 53] = 7.984000e+00;
   par[ 54] = 3.000000e+00;
   par[ 55] = 1.000000e+01;
   par[ 56] = 1.220600e+01;
   par[ 57] = 1.208000e+00;
   par[ 58] = 2.000000e-02;
   par[ 59] = 1.000000e+01;
   par[ 60] = 1.221000e+00;
   par[ 61] = 9.580000e-01;
   par[ 62] = 9.250000e-01;
   par[ 63] = 2.000000e+01;
   par[ 64] = 8.000000e-04;
   par[ 65] = 5.000000e+00;
   par[ 66] = 7.567000e-01;
   par[ 67] = 6.100000e-01;
   par[ 68] = 5.430000e-01;
   par[ 69] = 5.155800e+01;
   par[ 70] = 2.094500e+00;
   par[ 71] = 9.280000e+00;
   par[ 72] = 3.480270e+03;
   par[ 73] = 9.720000e-01;
   par[ 74] = 1.713710e+03;
   par[ 75] = 8.675140e+03;
   par[ 76] = 5.235000e+00;
   par[ 77] = 9.430000e-01;
   par[ 78] = 7.616400e+02;
   par[ 79] = 5.130000e+00;
   par[ 80] = 1.445000e+00;
   par[ 81] = 2.285000e+00;
   par[ 82] = 6.000000e+01;
   par[ 83] = 1.800000e+02;
   par[ 84] = 2.821100e+01;
   par[ 85] = 1.940700e+02;
   par[ 86] = 1.142500e+02;
   par[ 87] = 4.287000e+00;
   par[ 88] = 8.994300e+01;
   par[ 89] = 4.474700e+02;
   par[ 90] = 1.322402e+05;
   par[ 91] = 1.724500e+02;
   par[ 92] = 1.990000e-01;
   par[ 93] = 1.600000e+01;
   par[ 94] = 3.000000e+00;
   par[ 95] = 2.000000e+00;
   par[ 96] = 1.000000e+00;
   par[ 97] = 2.200000e+02;
   par[ 98] = 1.000000e+01;
   par[ 99] = 5.593000e-03;
   par[100] = 2.200000e+02;
   par[101] = 2.000000e+00;
   par[102] = 9.600000e+00;
   par[103] = 1.000000e+00;
   par[104] = 3.221800e+02;
   par[105] = 6.443500e+02;
   par[106] = 4.470000e-01;
   par[107] = 3.222000e+00;
   par[108] = 3.221800e+01;
   par[109] = 8.949000e-05;
   par[110] = 3.221800e+01;
   par[111] = 8.950000e-02;
   par[112] = 3.221800e+01;
   par[113] = 3.222000e+00;
   par[114] = 8.950000e-03;
 
 
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
   int ms = std::min(33,*ls);
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
      *ls = 33;
   }
 
//-----------------------------------------------------------------------
   int mp = std::min(115,*lp);
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
      *lp = 115;
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
   *nspe = 33;
   for (int j = 0; j < *nspe; ++j)
   {
      strncpy(ids[j],"\0",_len);
   }
   ///
   strncpy(ids[  0],"LH_pit\0",_len);
   strncpy(ids[  1],"LH_blood\0",_len);
   strncpy(ids[  2],"R_LH\0",_len);
   strncpy(ids[  3],"LH_R\0",_len);
   strncpy(ids[  4],"R_LH_des\0",_len);
   strncpy(ids[  5],"FSH_pit\0",_len);
   strncpy(ids[  6],"FSH_blood\0",_len);
   strncpy(ids[  7],"R_FSH\0",_len);
   strncpy(ids[  8],"FSH_R\0",_len);
   strncpy(ids[  9],"R_FSH_des\0",_len);
   strncpy(ids[ 10],"S_foll\0",_len);
   strncpy(ids[ 11],"AF1\0",_len);
   strncpy(ids[ 12],"AF2\0",_len);
   strncpy(ids[ 13],"AF3\0",_len);
   strncpy(ids[ 14],"AF4\0",_len);
   strncpy(ids[ 15],"PrF\0",_len);
   strncpy(ids[ 16],"OvF\0",_len);
   strncpy(ids[ 17],"Sc1\0",_len);
   strncpy(ids[ 18],"Sc2\0",_len);
   strncpy(ids[ 19],"Lut1\0",_len);
   strncpy(ids[ 20],"Lut2\0",_len);
   strncpy(ids[ 21],"Lut3\0",_len);
   strncpy(ids[ 22],"Lut4\0",_len);
   strncpy(ids[ 23],"E2\0",_len);
   strncpy(ids[ 24],"P4\0",_len);
   strncpy(ids[ 25],"InhA\0",_len);
   strncpy(ids[ 26],"InhB\0",_len);
   strncpy(ids[ 27],"InhA_tau\0",_len);
   strncpy(ids[ 28],"GnRH_G\0",_len);
   strncpy(ids[ 29],"R_GnRH_a\0",_len);
   strncpy(ids[ 30],"R_GnRH_i\0",_len);
   strncpy(ids[ 31],"GnRH_R_a\0",_len);
   strncpy(ids[ 32],"GnRH_R_i\0",_len);
   ///
}
//=======================================================================
void get_parameter_ids_ ( char (*idp)[MAXIDSTRLEN], int* npar, int _len)
{
   *npar = 115;
   for (int j = 0; j < *npar; ++j)
   {
      strncpy(idp[j],"\0",_len);
   }
   ///
   strncpy(idp[  0],"global_p_001_001\0",_len);
   strncpy(idp[  1],"global_p_001_002\0",_len);
   strncpy(idp[  2],"global_p_001_003\0",_len);
   strncpy(idp[  3],"global_p_001_004\0",_len);
   strncpy(idp[  4],"global_p_001_005\0",_len);
   strncpy(idp[  5],"global_p_001_006\0",_len);
   strncpy(idp[  6],"global_p_001_007\0",_len);
   strncpy(idp[  7],"global_p_001_008\0",_len);
   strncpy(idp[  8],"global_p_001_009\0",_len);
   strncpy(idp[  9],"global_p_001_010\0",_len);
   strncpy(idp[ 10],"global_p_002_001\0",_len);
   strncpy(idp[ 11],"global_p_002_002\0",_len);
   strncpy(idp[ 12],"global_p_002_003\0",_len);
   strncpy(idp[ 13],"global_p_003_001\0",_len);
   strncpy(idp[ 14],"global_p_004_001\0",_len);
   strncpy(idp[ 15],"global_p_006_001\0",_len);
   strncpy(idp[ 16],"global_p_006_002\0",_len);
   strncpy(idp[ 17],"global_p_006_003\0",_len);
   strncpy(idp[ 18],"global_p_006_004\0",_len);
   strncpy(idp[ 19],"global_p_006_005\0",_len);
   strncpy(idp[ 20],"global_p_006_006\0",_len);
   strncpy(idp[ 21],"global_p_006_007\0",_len);
   strncpy(idp[ 22],"global_p_006_008\0",_len);
   strncpy(idp[ 23],"global_p_006_009\0",_len);
   strncpy(idp[ 24],"global_p_006_010\0",_len);
   strncpy(idp[ 25],"global_p_006_011\0",_len);
   strncpy(idp[ 26],"global_p_007_001\0",_len);
   strncpy(idp[ 27],"global_p_007_002\0",_len);
   strncpy(idp[ 28],"global_p_008_001\0",_len);
   strncpy(idp[ 29],"global_p_009_001\0",_len);
   strncpy(idp[ 30],"global_p_011_001\0",_len);
   strncpy(idp[ 31],"global_p_011_002\0",_len);
   strncpy(idp[ 32],"global_p_011_003\0",_len);
   strncpy(idp[ 33],"global_p_011_004\0",_len);
   strncpy(idp[ 34],"global_p_011_005\0",_len);
   strncpy(idp[ 35],"global_p_011_006\0",_len);
   strncpy(idp[ 36],"global_p_011_007\0",_len);
   strncpy(idp[ 37],"global_p_012_001\0",_len);
   strncpy(idp[ 38],"global_p_012_002\0",_len);
   strncpy(idp[ 39],"global_p_012_003\0",_len);
   strncpy(idp[ 40],"global_p_012_004\0",_len);
   strncpy(idp[ 41],"global_p_013_001\0",_len);
   strncpy(idp[ 42],"global_p_013_002\0",_len);
   strncpy(idp[ 43],"global_p_013_003\0",_len);
   strncpy(idp[ 44],"global_p_014_001\0",_len);
   strncpy(idp[ 45],"global_p_014_002\0",_len);
   strncpy(idp[ 46],"global_p_014_003\0",_len);
   strncpy(idp[ 47],"global_p_014_004\0",_len);
   strncpy(idp[ 48],"global_p_015_001\0",_len);
   strncpy(idp[ 49],"global_p_015_002\0",_len);
   strncpy(idp[ 50],"global_p_015_003\0",_len);
   strncpy(idp[ 51],"global_p_016_001\0",_len);
   strncpy(idp[ 52],"global_p_016_002\0",_len);
   strncpy(idp[ 53],"global_p_017_001\0",_len);
   strncpy(idp[ 54],"global_p_017_002\0",_len);
   strncpy(idp[ 55],"global_p_017_003\0",_len);
   strncpy(idp[ 56],"global_p_017_004\0",_len);
   strncpy(idp[ 57],"global_p_018_001\0",_len);
   strncpy(idp[ 58],"global_p_018_002\0",_len);
   strncpy(idp[ 59],"global_p_018_003\0",_len);
   strncpy(idp[ 60],"global_p_018_004\0",_len);
   strncpy(idp[ 61],"global_p_019_001\0",_len);
   strncpy(idp[ 62],"global_p_020_001\0",_len);
   strncpy(idp[ 63],"global_p_020_002\0",_len);
   strncpy(idp[ 64],"global_p_020_003\0",_len);
   strncpy(idp[ 65],"global_p_020_004\0",_len);
   strncpy(idp[ 66],"global_p_021_001\0",_len);
   strncpy(idp[ 67],"global_p_022_001\0",_len);
   strncpy(idp[ 68],"global_p_023_001\0",_len);
   strncpy(idp[ 69],"global_p_024_001\0",_len);
   strncpy(idp[ 70],"global_p_024_002\0",_len);
   strncpy(idp[ 71],"global_p_024_003\0",_len);
   strncpy(idp[ 72],"global_p_024_004\0",_len);
   strncpy(idp[ 73],"global_p_024_005\0",_len);
   strncpy(idp[ 74],"global_p_024_006\0",_len);
   strncpy(idp[ 75],"global_p_024_007\0",_len);
   strncpy(idp[ 76],"global_p_024_008\0",_len);
   strncpy(idp[ 77],"global_p_025_001\0",_len);
   strncpy(idp[ 78],"global_p_025_002\0",_len);
   strncpy(idp[ 79],"global_p_025_003\0",_len);
   strncpy(idp[ 80],"global_p_026_001\0",_len);
   strncpy(idp[ 81],"global_p_026_002\0",_len);
   strncpy(idp[ 82],"global_p_026_003\0",_len);
   strncpy(idp[ 83],"global_p_026_004\0",_len);
   strncpy(idp[ 84],"global_p_026_005\0",_len);
   strncpy(idp[ 85],"global_p_026_006\0",_len);
   strncpy(idp[ 86],"global_p_026_007\0",_len);
   strncpy(idp[ 87],"global_p_026_008\0",_len);
   strncpy(idp[ 88],"global_p_027_001\0",_len);
   strncpy(idp[ 89],"global_p_027_002\0",_len);
   strncpy(idp[ 90],"global_p_027_003\0",_len);
   strncpy(idp[ 91],"global_p_027_004\0",_len);
   strncpy(idp[ 92],"global_p_028_001\0",_len);
   strncpy(idp[ 93],"global_p_029_001\0",_len);
   strncpy(idp[ 94],"global_p_029_002\0",_len);
   strncpy(idp[ 95],"global_p_029_003\0",_len);
   strncpy(idp[ 96],"global_p_029_004\0",_len);
   strncpy(idp[ 97],"global_p_029_005\0",_len);
   strncpy(idp[ 98],"global_p_029_006\0",_len);
   strncpy(idp[ 99],"global_p_030_001\0",_len);
   strncpy(idp[100],"global_p_030_002\0",_len);
   strncpy(idp[101],"global_p_030_003\0",_len);
   strncpy(idp[102],"global_p_030_004\0",_len);
   strncpy(idp[103],"global_p_030_005\0",_len);
   strncpy(idp[104],"global_p_031_001\0",_len);
   strncpy(idp[105],"global_p_031_002\0",_len);
   strncpy(idp[106],"global_p_031_003\0",_len);
   strncpy(idp[107],"global_p_032_001\0",_len);
   strncpy(idp[108],"global_p_032_002\0",_len);
   strncpy(idp[109],"global_p_033_001\0",_len);
   strncpy(idp[110],"global_p_033_002\0",_len);
   strncpy(idp[111],"global_p_033_003\0",_len);
   strncpy(idp[112],"global_p_034_001\0",_len);
   strncpy(idp[113],"global_p_034_002\0",_len);
   strncpy(idp[114],"global_p_035_001\0",_len);
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
   strncpy(idm[  0],"PAEON_V2\0",_len);
   strncpy(idm[  1],"Tue Jul 28 19:03:47 2015\0",_len);
   strncpy(idm[  2],"001438103027\0",_len);
   strncpy(idm[  3],"PAEON_V2.xml\0",_len);
   strncpy(idm[  4],"with vareq\0",_len);
   ///
}
//=======================================================================
};
//=======================================================================
//=======================================================================
template<typename T> 
      T fun0 ( T& x1, T& x2, T& x3 )
 
/*
      T&  x1;
      T&  x2;
      T&  x3;
 
*/
{
      T ret0;
      ret0 = pow(x1 / x2, x3) / (1.000000e+00 + pow(x1 / x2, x3));
 
      return ret0;
}
//=======================================================================
template<typename T> 
      T fun1 ( T& x1, T& x2, T& x3 )
 
/*
      T&  x1;
      T&  x2;
      T&  x3;
 
*/
{
      T ret1;
      ret1 = 1.000000e+00 / (1.000000e+00 + pow(x1 / x2, x3));
 
      return ret1;
}
 
 
//=======================================================================

