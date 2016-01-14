// For a description of this C code, see diagsqt.m.
#include "types.h"
#include "doublevectormex.h"
#include "singlematrixmex.h"
#include "misc.h"
#include "diagsq.h"
#include "mex.h"
#include "matrix.h"

// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  // Get the input matrix X.
  const SingleMatrix X = getSingleMatrix(prhs[0]);

  // Get the input vector a.
  const DoubleVector a = getDoubleVector(prhs[1]);
    
  // Initialize the output.
  DoubleVector y = createMatlabVector(X.nr,&plhs[0]);
  
  // Compute y = X.^2*a.
  diagsqt(X.elems,a.elems,y.elems,X.nr,X.nc);    
}
