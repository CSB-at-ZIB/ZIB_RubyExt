/*
c
c-----------------------------------------------------------------------
c
c            *************************************************
c            **                                             **
c            **                ydot_LIMEX                   **
c            **                                             **
c            *************************************************
c
c     ydot_LIMEX is an abstract interface
c
c
c-----------------------------------------------------------------------
c
c     Copyright (C) 2010 - 2014,
c     Zuse Institute Berlin (ZIB)
c     ALL RIGHTS RESERVED
c
c     Written by: automated transcription of 'sbml2fortran'
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

extern struct
{
    double start;
} sbmlvariables_;

/*
extern void get_species_ids_ ( char (*ids)[MAXIDSTRLEN], int* nspe, int _len = MAXIDSTRLEN );
extern void get_parameter_ids_ ( char (*idp)[MAXIDSTRLEN], int* npar, int _len = MAXIDSTRLEN );
extern void get_model_ids_ ( char (*idm)[MAXIDSTRLEN], int* nid, int _len = MAXIDSTRLEN );
*/
extern void get_species_ids_ ( char* ids, int* nspe, int _len );
extern void get_parameter_ids_ ( char* idp, int* npar, int _len );
extern void get_model_ids_ ( char* idm, int* nid, int _len );

extern void init_ode_ ( double* com, int* icom, int* lcom,
                        double* spe, int* ispe, int* lspe,
                        double* par, int* ipar, int* lpar );
extern void init_var_ ( int* nidx, int* pidx );

/* dummy routine here */
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

/* computing rhs f(t,y) and mass matrix B(t,y) */
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

/* computing rhs f(t,y), fy(t,y)S + fp(t,y), and related mass matrix B(t,y) */
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

