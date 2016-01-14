// For a description of this C code, see diagsq.m.
#include "types.h"
#include "doublevectormex.h"
#include "singlematrixmex.h"
#include "diagsq.h"
#include "mex.h"
#include "matrix.h"

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  // Get the input matrix X.
  const SingleMatrix X = getSingleMatrix(prhs[0]);

  // Get the input vector a.
  const DoubleVector a = getDoubleVector(prhs[1]);
  
  // Initialize the output.
  DoubleVector y = createMatlabVector(X.nc,&plhs[0]);
  
  // Compute y = (X.^2)'*a.
  diagsq(X.elems,a.elems,y.elems,X.nr,X.nc);
}
