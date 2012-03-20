// For a description of this C code, see diagsqt.m.
#include "types.h"
#include "doublevectormatlab.h"
#include "singlematrixmatlab.h"
#include "diagsqt.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// -----------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // GET INPUTS.
  // Get the input matrix X.
  const SingleMatrix X = getSingleMatrix(prhs[0]);

  // Get the dimensions of X.
  const Size n = X.nr;
  const Size p = X.nc;

  // Get the input vector a.
  const DoubleVector a = getDoubleVector(prhs[1]);
    
  // INITIALIZE THE OUTPUT.
  DoubleVector y = createMatlabVector(n,&plhs[0]);
  
  // COMPUTE y = X.^2*a.
  diagsqt(X.elems,a.elems,y.elems,X.nr,X.nc);    
}

