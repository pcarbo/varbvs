#include "singlematrix.h"

// Function definitions.
// -----------------------------------------------------------------
// Get information about a dense, single precision matrix from MATLAB.
singlematrix getsinglematrix (const mxArray* ptr) {
  singlematrix result;  // The return value.
  result.nr    = mxGetM(ptr);
  result.nc    = mxGetN(ptr);
  result.elems = (float*) mxGetPr(ptr);
  return result;
}

// Copy the jth column of matrix X into vector y.
void getcolumn (const singlematrix& X, doublevector& y, mwIndex j) {
  const mwSize m  = X.nr;           // Height of X.
  const float* xi = X.elems + m*j;  // Pointer to (i,j) entry of X.
  double*      yi = y.elems;        // Point to ith entry of y.

  // Repeat for each row.
  for (mwIndex i = 0; i < m; i++, xi++, yi++)
    *yi = (double) *xi;
}
