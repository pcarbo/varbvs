#include "doublematrixmatlab.h"

// Function definitions.
// -----------------------------------------------------------------
// Get a double-precision floating-point matrix from a MATLAB array. 
DoubleMatrix getDoubleMatrix (const mxArray* ptr) {
  DoubleMatrix result;
  result.nr    = mxGetM(ptr);
  result.nc    = mxGetN(ptr);
  result.elems = (double*) mxGetPr(ptr);
  return result;
}
