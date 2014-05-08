// For a description of this C code, see varbvsmixupdate.m.
#include "types.h"
// #include "vectorops.h"
#include "doublevectormatlab.h"
#include "doublematrixmatlab.h"
#include "varbvsmix.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// MEX-file gateway routine. Note that varbvsmixupdate.m checks the
// inputs, so we do not have to do it here.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // GET INPUTS.
  const DoubleMatrix X       = getDoubleMatrix(prhs[0]);
  const double       sigma   = *mxGetPr(prhs[1]);
  const double       sa1     = *mxGetPr(prhs[2]);
  const double       sa2     = *mxGetPr(prhs[3]);
  const DoubleVector logodds = getDoubleVector(prhs[4]);
  const DoubleVector xy      = getDoubleVector(prhs[5]);
  const DoubleVector d       = getDoubleVector(prhs[6]);
  const DoubleVector alpha0  = getDoubleVector(prhs[7]);
  const DoubleVector mu10    = getDoubleVector(prhs[8]);
  const DoubleVector mu20    = getDoubleVector(prhs[9]);
  const DoubleVector Xr0     = getDoubleVector(prhs[10]);
  const DoubleVector I       = getDoubleVector(prhs[11]);

  // Get the number of samples (n), the number of variables (p), and
  // the number of coordinate ascent updates (m).
  const Size n = X.nr;
  const Size p = X.nc;
  const Size m = I.n;

  // INITIALIZE OUTPUTS.
  DoubleVector alpha = createMatlabVector(p,&plhs[0]);
  DoubleVector mu1   = createMatlabVector(p,&plhs[1]);
  DoubleVector mu2   = createMatlabVector(p,&plhs[2]);
  DoubleVector Xr    = createMatlabVector(n,&plhs[3]);

  copyDoubleVector(alpha0,alpha);
  copyDoubleVector(mu10,mu1);
  copyDoubleVector(mu20,mu2);
  copyDoubleVector(Xr0,Xr);

  // This is storage for a column of matrix X.
  double* x = malloc(sizeof(double)*n);

  // RUN COORDINATE ASCENT UPDATES.
  // Repeat for each coordinate ascent update.
  for (Index j = 0; j < m; j++) {
    Index k = (Index) I.elems[j];

    // Copy the kth column of matrix X.
    copyDoubleColumn(X.elems,x,k,n);

    // Perform the update.
    varbvsmixupdate(x,xy.elems[k],d.elems[k],sigma,sa1,sa2,logodds.elems[k],
		    alpha.elems+k,mu1.elems+k,mu2.elems+k,Xr.elems,n);
  }

  // Free the dynamically allocated memory.
  free(x);
}
