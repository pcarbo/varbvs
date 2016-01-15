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
