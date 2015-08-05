// Copyright (C) 2010 - 2015
// ZIB - Zuse Institute Berlin, Germany
//
// first added : 2015-06-30 td
// last changed:
//
#ifndef __YDOT_LIMEX_CPP_H
#define __YDOT_LIMEX_CPP_H

/*
c
c-----------------------------------------------------------------------
c
c            *************************************************
c            **                                             **
c            **                ydot_LIMEXcpp                **
c            **                                             **
c            *************************************************
c
c     ydot_LIMEXcpp is an abstract interface
c
c
c-----------------------------------------------------------------------
c
c     Copyright (C) 2010 - 2015,
c     Zuse Institute Berlin (ZIB)
c     ALL RIGHTS RESERVED
c
c     Written by: automated transcription of 'sbml2adolc'
c
c     Zuse Institute Berlin (ZIB)
c     Takustrasse 7
c     D-14195 Berlin-Dahlem, Germany
c
c     Phone : +49-30-84185-0
c     Fax   : +49-30-84185-125
c     URL   : http://www.zib.de
c
*/
#define MAXIDSTRLEN 64


extern "C" 
{
  extern struct
  {
      double start;
  } sbmlvariables_;



extern void get_species_ids_ ( char (*ids)[MAXIDSTRLEN], int* nspe, 
                               int _len = MAXIDSTRLEN );
extern void get_parameter_ids_ ( char (*idp)[MAXIDSTRLEN], int* npar, 
                                 int _len = MAXIDSTRLEN );
extern void get_model_ids_ ( char (*idm)[MAXIDSTRLEN], int* nid, 
                             int _len = MAXIDSTRLEN );


extern void init_ode_ ( double* com, int* icom, int* lcom,
                        double* spe, int* ispe, int* lspe,
                        double* par, int* ipar, int* lpar );

extern void init_var_ ( int* nidx, int* pidx );

extern void set_model_number_directions_ ();

// dummy routine here
/*
extern void ydot_jacobian_ (
                    int*        n,
                    double*     t,
                    double*     y,
                    double*     dy,
                    double*     J,
                    int*        ldJ,
                    int*        ml,
                    int*        mu,
                    int*        full_or_band,
                    int*        info
             );
*/

// computing rhs f(t,y) and mass matrix B(t,y)
extern void ydot_limex_ (
                    int*        n,
                    int*        nz,
                    double*     t,
                    double*     y,
                    double*     dy,
                    double*     B,
                    int*        ir,
                    int*        ic,
                    int*        info
             );


// computing rhs f(t,y), fy(t,y)S + fp(t,y), and the related mass matrix B(t,y)
extern void ydot_slimex_ (
                    int*        n,
                    int*        nz,
                    double*     t,
                    double*     y,
                    double*     dy,
                    double*     B,
                    int*        ir,
                    int*        ic,
                    int*        info
             );
}

#endif // __YDOT_LIMEX_CPP_H
