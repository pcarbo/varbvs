#include "doublematrixmex.h"

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// Get a double precision matrix from a MATLAB array. 
DoubleMatrix getDoubleMatrix (const mxArray* ptr) {
  DoubleMatrix result;
  result.nr    = mxGetM(ptr);
  result.nc    = mxGetN(ptr);
  result.elems = (double*) mxGetPr(ptr);
  return result;
}

// Get a pointer to column j of matrix X.
const double* getDoubleColumn (const double* X, Index j, Size n) {
  return X + n*j;
}

// Copy column j of matrix X.
void copyDoubleColumn (const double* X, double* y, Index j, Size n) {
  const double* xij = getDoubleColumn(X,j,n);
  for (Index i = 0; i < n; i++, xij++, y++)
    *y = (double) *xij;
}

