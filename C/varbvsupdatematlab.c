// For a description of this C code, see varbvsupdate.m.
#include "doublevectormatlab.h"
#include "singlematrixmatlab.h"
#include "vectorops.h"
#include "varbvs.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// -----------------------------------------------------------------
// MEX-file gateway routine. Note that varbvsupdate.m checks the
// inputs, so we do not have to do it here.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // GET INPUTS.
  // Get the input matrix X.
  const SingleMatrix X = getSingleMatrix(prhs[0]);

  // Get input scalars sigma and sa.
  const double sigma = *mxGetPr(prhs[1]);
  const double sa    = *mxGetPr(prhs[2]);

  // Get input vector logodds.
  const DoubleVector logodds = getDoubleVector(prhs[3]);

  // Get input vectors xy and d.
  const DoubleVector xy = getDoubleVector(prhs[4]);
  const DoubleVector d  = getDoubleVector(prhs[5]);

  // Get input vectors alpha0 and mu0.
  const DoubleVector alpha0 = getDoubleVector(prhs[6]);
  const DoubleVector mu0    = getDoubleVector(prhs[7]);

  // Get the input vector Xr0.
  const DoubleVector Xr0 = getDoubleVector(prhs[8]);

  // Get input vector I.
  const DoubleVector I = getDoubleVector(prhs[9]);

  // Get the number of samples (n) and the number of variables (p).
  const Size n = X.nr;
  const Size p = X.nc;

  // Get the number of coordinate ascent updates.
  Size m = I.n;

  // INITIALIZE OUTPUTS.
  DoubleVector alpha = createMatlabVector(p,&plhs[0]);
  DoubleVector mu    = createMatlabVector(p,&plhs[1]);
  DoubleVector Xr    = createMatlabVector(n,&plhs[2]);

  copyDoubleVector(alpha0,alpha);
  copyDoubleVector(mu0,mu);
  copyDoubleVector(Xr0,Xr);

  // This is storage for a column of matrix X.
  double* x = malloc(sizeof(double)*n);

  // RUN COORDINATE ASCENT UPDATES.
  // Repeat for each coordinate ascent update.
  for (Index j = 0; j < m; j++) {

    // We need to subtract 1 from the index here because MATLAB arrays
    // start at one, but C arrays start at zero.
    Index k = (Index) I.elems[j] - 1;

    // Get the kth column of matrix X.
    getColumn(X.elems,x,k,X.nr);

    // Perform the update.
    varbvsupdate(x,xy.elems[k],d.elems[k],sigma,sa,logodds.elems[k],
		 alpha.elems+k,mu.elems+k,Xr.elems,n);
  }

  // Free the dynamically allocated memory.
  free(x);
}
