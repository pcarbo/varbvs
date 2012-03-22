// For a description of this C code, see diagsq.m.
#include "types.h"
#include "doublevectormatlab.h"
#include "singlematrixmatlab.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Function declarations.
// -----------------------------------------------------------------
// Compute (X.^2)'*a and store the result in vector y. X is an m x n
// matrix in which the entries in each column are stored consecutively
// in memory.
void diagsq (const MatrixElem* X, const double* a, double* y, 
	     Size m, Size n);

// Function definitions.
// -----------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // Get the input matrix X.
  const SingleMatrix X = getSingleMatrix(prhs[0]);

  // Get the input vector a.
  const DoubleVector a = getDoubleVector(prhs[1]);
  
  // Initialize the output.
  DoubleVector y = createMatlabVector(X.nc,&plhs[0]);
  
  // Compute y = (X.^2)'*a.
  diagsq(X.elems,a.elems,y.elems,X.nr,X.nc);
}

// Compute (X.^2)'*a and store the result in vector y.
void diagsq (const MatrixElem* X, const double* a, double* y, 
	     Size m, Size n) {
  double t;

  // Repeat for each column of X.
  for (Index j = 0; j < n; j++, y++) {
    
    // Initialize the jth entry of the result.
    *y = 0;

    // Repeat for each row of X.
    const double* ai = a;
    for (Index i = 0; i < m; i++, X++, ai++) {

      // Add X(i,j)^2 * a(i) to the jth entry of y.
      t   = (double) *X;
      *y += t * t * (*ai);
    }
  }
}
