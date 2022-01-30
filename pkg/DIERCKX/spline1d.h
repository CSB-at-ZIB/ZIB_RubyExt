/* 
*/

typedef struct
  { 
	double	*t;
	int	n;	/* n = length(t) */
	double	*coeffs;
	int	k;	/* order k: 1 <= k <= 5 */
	int	bc;	/* boundary cond: 0="extrapolate", 1="zero", 2="error", 3="nearest" */
	double	residual;
	double	*wrk;
	int	lwrk;
	int	*iwrk;
	int	ierror;
  }
spline1d_state_t;


/* DIERCKX API for Spline1D */

extern void curfit_(
	int*	iopt,
	int*	m,
	double*	xin,
	double*	yin,
	double*	win,
	double*	xb,
	double*	xe,
	int*	k,
	double* s,
	int*	nest,
	int*	n,
	double*	t,
	double*	coeffs,
	double*	residual,
	double*	wrk,
	int*	lwk,
	int*	iwrk,
	int*	ierror);


extern void percur_(
	int*	iopt,
	int*	m,
	double*	xin,
	double*	yin,
	double*	win,
	int*	k,
	double* s,
	int*	nest,
	int*	n,
	double*	t,
	double*	coeffs,
	double*	residual,
	double*	wrk,
	int*	lwk,
	int*	iwrk,
	int*	ierror);


extern void splev_(
	double*	t,
	int*	n,
	double* coeffs,
	int*	k,
	double* xin,
	double*	yout,
	int*	m,
	int*	bc,
	int*	ierror);


extern void splder_(
	double*	t,
	int*	n,
	double* coeffs,
	int*	k,
	int*	nu,
	double* xin,
	double*	yout,
	int*	m,
	int*	bc,
	double*	wrk,
	int*	ierror);


extern double* splint_(
	double*	t,
	int*	n,
	double* coeffs,
	int*	k,
	double* a,
	double*	b,
	double*	wrk);


extern void sproot_(
	double*	t,
	int*	n,
	double* coeffs,
	double* zeros,
	int*	mest,
	int*	m,
	int*	ierror);

