// For a description of this C++ code, see var1fast.m
#include "common.h"
#include "doublevector.h"
#include "singlematrix.h"

// These include files have a bunch of definitions to interface C++
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Function declarations.
// -----------------------------------------------------------------
// Compute the sample mean.
double computemean (const doublevector& x);

// Compute the sample variance.
double computevariance (const doublevector& x);

// -----------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {
  const mxArray* ptr = prhs[0];

  // Make sure we have the correct number of inputs and outputs.
  if (nrhs != 1)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 1)
    mexErrMsgTxt("Incorrect number of output arguments");

  // Get the input matrix X.
  if (mxGetClassID(ptr) != mxSINGLE_CLASS)
    mexErrMsgTxt("Input argument X must be SINGLE");
  const singlematrix X = getsinglematrix(ptr);

  // Get the number of samples (n) and the number of variables (p).
  const int n = X.nr;
  const int p = X.nc;

  // Initialize the output.
  doublevector y = creatematlabvector(p,plhs[0]);

  // This is storage for columns of X.
  doublevector x = newdoublevector(n);

  // Repeat for each column of X.
  for (int k = 0; k < p; k++) {

    // Get the kth column of X.
    getcolumn(X,x,k);

    // Compute the sample variance.
    y.elems[k] = computevariance(x);
  }

  // Free dynamically allocated memory.
  free(x);
}

// Compute the sample mean.
double computemean (const doublevector& x) {
  return sum(x) / x.n;
}

// Compute the sample variance.
double computevariance (const doublevector& x) {
  double y = 0;  // The return value.
  double t;      // Intermediate result.

  // Compute the sample mean.
  double mu = computemean(x);

  // Repeat for each entry of x.
  for (int i = 0; i < x.n; i++) {
    t  = x.elems[i] - mu;
    y += t*t;
  }

  // Divide by the number of entries.
  y /= x.n;

  return y;
}
