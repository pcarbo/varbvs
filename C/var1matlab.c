// For a description of this C code, see var1.m.
#include "types.h"
#include "vectorops.h"
#include "doublevectormatlab.h"
#include "singlematrixmatlab.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Function declarations.
// -----------------------------------------------------------------
// Compute the sample variance.
double computeVar (const double* x, int n);

// -----------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // Get the input matrix X.
  const SingleMatrix X = getSingleMatrix(prhs[0]);

  // Get the number of samples (n) and the number of variables (p).
  const int n = X.nr;
  const int p = X.nc;

  // Initialize the output.
  DoubleVector y = createMatlabVector(p,&plhs[0]);

  // This is storage for a column of matrix X.
  double* xk = malloc(sizeof(double)*n);

  // Repeat for each column of X.
  for (int k = 0; k < p; k++) {

    // Get the kth column of X.
    copyColumn(X.elems,xk,k,n);

    // Compute the sample variance.
    y.elems[k] = computeVar(xk,n);
  }

  // Free the dynamically allocated memory.
  free(xk);
}

// Compute the sample variance.
double computeVar (const double* x, int n) {
  double y = 0;  // The return value.
  double t;      // Intermediate result.

  // Compute the sample mean.
  double mu = sum(x,n)/n;

  // Repeat for each entry of x.
  for (int i = 0; i < n; i++) {
    t  = x[i] - mu;
    y += t*t;
  }

  // Divide by the number of entries.
  y /= n;

  return y;
}
