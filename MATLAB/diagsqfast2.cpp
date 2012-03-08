// For a description of this C++ code, see diagsqfast2.m.
#include "mex.h"
#include "matrix.h"
#include "doublevector.h"
#include "singlematrix.h"

// Macros.
// ---------------------------------------------------------------
// Return the square of x.
#define square(x) ((x)*(x))

// Function declarations.
// -----------------------------------------------------------------
// Compute X.^2*A and store the result in Y.
void diagsq2 (const singlematrix& X, const doublevector& a, doublevector& y);

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
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument A must be a double precision vector \
of length P");
  const doublevector a = getdoublevector(ptr);
    
  // (2.) INITIALIZE OUTPUT.
  doublevector y = creatematlabvector(n,plhs[0]);
  
  // (3.) COMPUTE Y = X.^2*A.
  diagsq2(X,a,y);    
}

// Compute X.^2*A and store the result in Y.
void diagsq2 (const singlematrix& X, const doublevector& a, doublevector& y) {
  mwSize  m = X.nr;  // Size of left-hand vector.
  mwSize  n = X.nc;  // Size of right-hand vector.
  mwIndex i;         // Row of X.
  mwIndex j;         // Column of X.
  double* yi;        // Pointer to ith element of vector y.
  double* aj;        // Pointer to jth element of vector a.
  float*  xij;       // Pointer to (i,j) element of matrix X.
  double  t;         // An intermediate result.

  // Zero the entries of the result vector.
  setvector(y,0);

  // Repeat for each column of X.
  for (j = 0, xij = X.elems, aj = a.elems; j < n; j++, aj++)

    // Repeat for each row of X.
    for (i = 0, yi = y.elems; i < m; i++, xij++, yi++) {
      
      // Add X(i,j)^2 * a(j) to the ith entry of y.
      t    = (double) *xij;
      *yi += square(t) * (*aj);
    }
}
