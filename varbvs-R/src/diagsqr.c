// Part of the varbvs package, https://github.com/pcarbo/varbvs
//
// Copyright (C) 2012-2017, Peter Carbonetto
//
// This program is free software: you can redistribute it under the
// terms of the GNU General Public License; either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANY; without even the implied warranty of
// MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
#include <R.h>
#include <Rinternals.h>
#include "diagsq.h"

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// This function is used by diagsq in misc.R. This function is called
// in R using the .Call interface.
SEXP diagsq_Call (SEXP Xp, SEXP ap, SEXP yp) {

  // Get the inputs and outputs.
  double* X = REAL(Xp);  // Input X.
  double* a = REAL(ap);  // Input a.
  double* y = REAL(yp);  // Output y.

  // Get the number of rows (m) and columns (n) of X.
  SEXP     d = getAttrib(Xp,R_DimSymbol);
  R_xlen_t m = INTEGER(d)[0];
  R_xlen_t n = INTEGER(d)[1];

  // Compute y = (X.^2)'*a.
  diagsq(X,a,y,m,n);

  return R_NilValue;
}

// -----------------------------------------------------------------
// This function is used by diagsqt in misc.R. This function is called
// in R using the .Call interface.
SEXP diagsqt_Call (SEXP Xp, SEXP ap, SEXP yp) {

  // Get the inputs and outputs.
  double* X = REAL(Xp);  // Input X.
  double* a = REAL(ap);  // Input a.
  double* y = REAL(yp);  // Output y.

  // Get the number of rows (m) and columns (n) of X.
  SEXP     d = getAttrib(Xp,R_DimSymbol);
  R_xlen_t m = INTEGER(d)[0];
  R_xlen_t n = INTEGER(d)[1];

  // Compute y = (X.^2)'*a.
  diagsqt(X,a,y,m,n);

  return R_NilValue;
}
