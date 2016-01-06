#include "singlematrixmex.h"

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// Get a single precision matrix from a MATLAB array. 
SingleMatrix getSingleMatrix (const mxArray* ptr) {
  SingleMatrix result;
  result.nr = mxGetM(ptr);
  result.nc = mxGetN(ptr);
  result.elems = (MatrixElem*) mxGetPr(ptr);
  return result;
}
