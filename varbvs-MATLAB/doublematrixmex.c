#include "doublematrixmex.h"
#include "misc.h"

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

// Create an m x n matrix in a MATLAB array.
DoubleMatrix createMatlabDoubleMatrix (Size m, Size n, mxArray** ptr) {
  *ptr = mxCreateDoubleMatrix(m,n,mxREAL);
  return getDoubleMatrix(*ptr);
}

// Copy all the entries of one matrix to another matrix.
void copyDoubleMatrix (const DoubleMatrix source, DoubleMatrix dest) {
  copy(source.elems,dest.elems,dest.nr*dest.nc);
}

// Get a pointer to column j of matrix X.
const double* getConstDoubleColumn (const double* X, Index j, Size n) {
  return X + n*j;
}

// Get a pointer to column j of matrix X.
double* getDoubleColumn (double* X, Index j, Size n) {
  return X + n*j;
}

// Copy column j of matrix X.
void copyDoubleColumn (const double* X, double* y, Index j, Size n) {
  const double* xij = getConstDoubleColumn(X,j,n);
  for (Index i = 0; i < n; i++, xij++, y++)
    *y = (double) *xij;
}

