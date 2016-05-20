/*
 ============================================================================
 Name        : utility.c
 Author      : Evan Ray
 Version     :
 Copyright   : 
 Description : utility functions for numeric calculations and interacting with R
 ============================================================================
 */

#include <stdio.h>
#include <stdlib.h>

#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/Applic.h>
#include <R_ext/Lapack.h>

#include <Rdefines.h>
#include <R_ext/Rdynload.h>

#include <math.h>
#include <float.h>

#include "utility.h"

// Return the largest double that C can deal with
SEXP get_dbl_max() {
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	*(REAL(retval)) = DBL_MAX;
	
	UNPROTECT(1);
	return retval;
}

// Return the element with name matching str from the list provided.
// Based on http://cran.r-project.org/doc/manuals/R-exts.html#Handling-lists
//
// list is an R list object
// str is a character vector with the name of the list component to get
//
// Return a SEXP with the specified list component
SEXP getListElement(SEXP list, const char *str) {
	SEXP elmt = R_NilValue, names = getAttrib(list, R_NamesSymbol);

	for (R_len_t i = 0; i < length(list); i++)
		if(strcmp(CHAR(STRING_ELT(names, i)), str) == 0) {
			elmt = VECTOR_ELT(list, i);
			break;
		}

	return elmt;
}


// An interface to R's C function logspace_add
// Computes log(exp(log_x) + exp(log_y))
// Unlike R's implementation, we return
// -INFINITY if both log_x and log_y are -INFINITY
//
// log_x, log_y are length 1 numeric R vectors
//
// Return a SEXP with the length 1 numeric computed value
SEXP logspace_add_C(SEXP log_x, SEXP log_y) {
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	double lx = *(REAL(log_x)), ly = *(REAL(log_y));
	
	*(REAL(retval)) = logspace_add_safe(lx, ly);
	
	UNPROTECT(1);
	return retval;
}


// An interface to R's C function logspace_add
// Computes log(exp(log_x) + exp(log_y))
// Unlike R's implementation, we return
// -INFINITY if both log_x and log_y are -INFINITY
// 
// log_x, log_y are doubles
//
// Return a double with the computed value
double logspace_add_safe(double log_x, double log_y) {
	if(log_x == -INFINITY && log_y == -INFINITY) {
		return(-INFINITY);
	} else {
		return(logspace_add(log_x, log_y));
	}
}

// An interface to R's C function logspace_sub
// Computes log(exp(log_x) - exp(log_y))
//
// log_x, log_y are length 1 numeric R vectors
//
// Return a SEXP with the length 1 numeric computed value
SEXP logspace_sub_C(SEXP log_x, SEXP log_y) {
	SEXP retval;
	retval = PROTECT(allocVector(REALSXP, 1));
	
	double lx = *(REAL(log_x)), ly = *(REAL(log_y));

	*(REAL(retval)) = logspace_sub(lx, ly);
	
	UNPROTECT(1);
	return retval;
}

// Given a matrix X, compute the row-sums of X in log space
//
// Xp is an R matrix object with N_rowp rows and N_colp columns
//
// Return a SEXP with a numeric R vector of length N_rowp
SEXP logspace_sum_matrix_rows_C(SEXP Xp, SEXP N_rowp, SEXP N_colp) {
	int i, j, n_row = *INTEGER(N_rowp), n_col = *INTEGER(N_colp);
	SEXP retval = PROTECT(allocVector(REALSXP, n_row));
	double *dblptr = REAL(retval), *X = REAL(Xp);
	
	for(i = 0; i < n_row; i++) {
		*(dblptr + i) = *(X + i);
	}
	
	for(j = 1; j < n_col; j++) {
		for(i = 0; i < n_row; i++) {
			if(!(*(dblptr + i) == R_NegInf && *(X + i + j*n_row) == R_NegInf))
				*(dblptr + i) = logspace_add_safe(*(dblptr + i), *(X + i + j*n_row));
		}
	}
	
	UNPROTECT(1);
	return retval;
}

// Given a matrix X, compute the logspace-difference between the first and
// second columns: log(exp(X[i, 1]) - exp(X[i, 2])) for i = 1, ..., N_rowp
//
// Xp is an R matrix object with N_rowp rows and 2 columns
//
// Return a SEXP with a numeric R vector of length N_rowp
SEXP logspace_sub_matrix_rows_C(SEXP Xp, SEXP N_rowp) {
	int i, n_row = *INTEGER(N_rowp);
	SEXP retval = PROTECT(allocVector(REALSXP, n_row));
	double *dblptr = REAL(retval), *X = REAL(Xp);
	
	for(i = 0; i < n_row; i++) {
		*(dblptr + i) = logspace_sub(*(X + i), *(X + i + n_row));
	}

	UNPROTECT(1);
	return retval;
}

// Given a matrix Xp with observation vectors in the rows,
// determine whether each element of each observation vector is
// within the corresponding bounds specified by lowerp and upperp.
// 
// Xp is an r by c R numeric matrix object
// lowerp is an R numeric vector of length c specifying lower bounds for the
//   truncation region
// upperp is an R numeric vector of length c specifying upper bounds for the
//   truncation region
//
// Return a SEXP with an R integer vector of length r: element i is 1 if
//   lowerp[j] <= X[i, j] <= upperp[j] for all j = 1, ..., c AND
//   X[i, j] != -INFINITY AND X[i, j] != INFINITY,
//   and element i is 0 otherwise
SEXP in_truncation_support_C(SEXP Xp, SEXP lowerp, SEXP upperp) {
	int i, j,
		n_row = *(INTEGER(getAttrib(Xp, R_DimSymbol))),
		n_col = *(INTEGER(getAttrib(Xp, R_DimSymbol)) + 1);

	SEXP retval = PROTECT(allocVector(INTSXP, n_row));
	int *intptr = INTEGER(retval);
	
	double *X = REAL(Xp),
		*lower = REAL(lowerp),
		*upper = REAL(upperp);
	
	for(i = 0; i < n_row; i++) {
		*(intptr + i) = 1;
		for(j = 0; j < n_col; j++) {
			if(*(X + i + j*n_row) < *(lower + j) ||
				*(X + i + j*n_row) > *(upper + j) ||
				*(X + i + j*n_row) == R_NegInf ||
				*(X + i + j*n_row) == R_PosInf) {
				*(intptr + i) = 0;
				break;
			}
		}
	}
	
	UNPROTECT(1);
	return retval;
}


// Stuff below here lets R access the specified C functions
R_CallMethodDef callMethods[] =
{
    {"logspace_add_C", (DL_FUNC)&logspace_sub_C, 2},
    {"logspace_sum_matrix_rows_C", (DL_FUNC)&logspace_sum_matrix_rows_C, 3},
    {"logspace_sub_C", (DL_FUNC)&logspace_sub_C, 2},
    {"logspace_sub_matrix_rows_C", (DL_FUNC)&logspace_sub_matrix_rows_C, 2},
    {"in_truncation_support_C", (DL_FUNC)&in_truncation_support_C, 3},
    {NULL,NULL, 0}
};

void R_init_pdtmvn(DllInfo *dll)
{
    R_registerRoutines(dll,NULL,callMethods,NULL,NULL);
}

