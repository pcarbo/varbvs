#include "doublevectormex.h"
#include "misc.h"

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// Get a double precision vector from a MATLAB array.
DoubleVector getDoubleVector (const mxArray* ptr) {
  DoubleVector result;
  result.n = mxGetNumberOfElements(ptr);
  result.elems = mxGetPr(ptr);
  return result;
}

// Create a column vector in MATLAB array.
DoubleVector createMatlabVector (Size n, mxArray** ptr) {
  *ptr = mxCreateDoubleMatrix(n,1,mxREAL);
  return getDoubleVector(*ptr);
}

// Copy entries of one vector to another vector.
void copyDoubleVector (const DoubleVector source, DoubleVector dest) {
  copy(source.elems,dest.elems,dest.n);
}
