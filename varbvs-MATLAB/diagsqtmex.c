// For a description of this C code, see diagsqt.m.
#include "types.h"
#include "doublevectormex.h"
#include "singlematrixmex.h"
#include "vectorops.h"
#include "mex.h"
#include "matrix.h"

// FUNCTION DECLARATIONS
// ---------------------------------------------------------------------
// Compute X.^2*a and store the result in vector y. X is an m x n
// matrix in which the entries in each column are stored consecutively
// in memory.
void diagsqt (const MatrixElem* X, const double* a, double* y, Size m, Size n);

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

// ---------------------------------------------------------------------
// Compute X.^2*a and store the result in vector y.
void diagsqt (const MatrixElem* X, const double* a, double* y, Size m,
	      Size n) {
  double t;  // An intermediate result.

  // Zero the entries of the result vector.
  setVector(y,m,0);

  // Repeat for each column of X.
  for (Index j = 0; j < n; j++, a++) {

    // Repeat for each row of X.
    double* yi = y;
    for (Index i = 0; i < m; i++, X++, yi++) {
      
      // Add X(i,j)^2 * a(j) to the ith entry of y.
      t    = (double) *X;
      *yi += t * t * (*a);
    }
  }
}
