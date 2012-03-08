// For a description of this C++ code, see diagsqfast.m.
#include "doublevector.h"
#include "singlematrix.h"

// These include files have a bunch of definitions to interface C++
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Macros.
// ---------------------------------------------------------------
// Return the square of x.
#define square(x) ((x)*(x))

// Function declarations.
// -----------------------------------------------------------------
// Compute (X.^2)'*A and store the result in Y.
void diagsq (const singlematrix& X, const doublevector& a, doublevector& y);

// -----------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {
  const mxArray* ptr = prhs[0];

  // Make sure we have the correct number of inputs and outputs.
  if (nrhs != 2)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 1)
    mexErrMsgTxt("Incorrect number of output arguments");

  // (1.) GET INPUTS.
  // Get the input matrix X.
  if (mxGetClassID(ptr) != mxSINGLE_CLASS)
    mexErrMsgTxt("Input argument X must be SINGLE");
  const singlematrix X = getsinglematrix(ptr);

  // Get the size of the data set (n) and the number of SNPs (p).
  const mwSize n = X.nr;
  const mwSize p = X.nc;

  // Get the input vector a.
  ptr = prhs[1];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != n)
    mexErrMsgTxt("Input argument A must be a double precision vector \
of length N");
  const doublevector a = getdoublevector(ptr);
  
  // (2.) INITIALIZE OUTPUT.
  doublevector y = creatematlabvector(p,plhs[0]);
  
  // (3.) COMPUTE Y = (X.^2)'*a.
  diagsq(X,a,y);
}

// Compute (X.^2)'*A and store the result in Y.
void diagsq (const singlematrix& X, const doublevector& a, doublevector& y) {
  mwSize  m = X.nr;  // Size of left-hand vector.
  mwSize  n = X.nc;  // Size of right-hand vector.
  mwIndex i;         // Row of X.
  mwIndex j;         // Column of X.
  double* ai;        // Pointer to ith element of vector a.
  double* yj;        // Pointer to jth element of vector y.
  float*  xij;       // Pointer to (i,j) element of matrix X.
  double  t;         // An intermediate result.
  double  z;         // An intermediate result.

  // Repeat for each column of X.
  for (j = 0, xij = X.elems, yj = y.elems; j < n; j++, yj++) {
    
    // Initialize the jth entry of the result.
    z = 0;

    // Repeat for each row of X.
    for (i = 0, ai = a.elems; i < m; i++, xij++, ai++) {
      
      // Add X(i,j)^2 * a(i) to the jth entry of y.
      t  = (double) *xij;
      z += square(t) * (*ai);
    }

    // Store the result.
    *yj = z;
  }
}
