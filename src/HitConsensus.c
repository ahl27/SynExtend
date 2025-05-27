/****************************************************************************
 * calculate the offset between hit bounds and feature bounds in the linked 
 * pairs
 * author: nicholas cooley
 ****************************************************************************/

/*
 * Rdefines.h is needed for the SEXP typedef, for the error(), INTEGER(),
 * GET_DIM(), LOGICAL(), NEW_INTEGER(), PROTECT() and UNPROTECT() macros,
 * and for the NA_INTEGER constant symbol.
 */
#include <Rdefines.h>
#include "SEutils.h"

/*
 * R_ext/Rdynload.h is needed for the R_CallMethodDef typedef and the
 * R_registerRoutines() prototype.
 */
#include <R_ext/Rdynload.h>

/* for R_CheckUserInterrupt */
#include <R_ext/Utils.h>

// compute difference between hits
SEXP HitConsensus(SEXP f1lb, SEXP f1rb, SEXP f2lb, SEXP f2rb, SEXP s1, SEXP s2)
{
	int m1;
	int l01 = length(f1lb);
	
	double *feat1lb = REAL(f1lb);
	double *feat1rb = REAL(f1rb);
	double *feat2lb = REAL(f2lb);
	double *feat2rb = REAL(f2rb);
	int *strand1 = LOGICAL(s1);
	int *strand2 = LOGICAL(s2);
	
	// instantiate our output S expression
	SEXP res;
	PROTECT(res = allocVector(REALSXP, l01));
	// use a pointer to fill the vector with results
	double *resptr = REAL(res);
	
	for (m1 = 0; m1 < l01; m1++) {
		double temp1, temp2;
		if (strand1[m1] == strand2[m1]) {
			temp1 = feat1lb[m1] - feat2lb[m1];
			temp2 = feat1rb[m1] - feat2rb[m1];
		} else {
			temp1 = feat1lb[m1] - feat2rb[m1];
			temp2 = feat1rb[m1] - feat2lb[m1];
		}
		if (temp1 < 0)
			temp1 *= -1;
		if (temp2 < 0)
			temp2 *= -1;
		resptr[m1] = (temp1 + temp2)/2;
	}
	
	// unprotection occurs on a first in first out basis?
	UNPROTECT(1);
	
	return res;
}
