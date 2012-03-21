// For a description of this C code, see varbvsbinupdate.m.
#include "types.h"
#include "vectorops.h"
#include "doublevectormatlab.h"
#include "singlematrixmatlab.h"
#include "varbvsbin.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// MEX-file gateway routine. Note that varbvsbinupdate.m checks the
// inputs, so we do not have to do it here.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // GET INPUTS.
  const SingleMatrix X       = getSingleMatrix(prhs[0]);
  const double       sa      = *mxGetPr(prhs[1]);
  const DoubleVector logodds = getDoubleVector(prhs[2]);
  const DoubleVector u       = getDoubleVector(mxGetField(prhs[3],0,"u"));
  const DoubleVector d       = getDoubleVector(mxGetField(prhs[3],0,"d"));
  const DoubleVector xy      = getDoubleVector(mxGetField(prhs[3],0,"xy"));
  const DoubleVector xu      = getDoubleVector(mxGetField(prhs[3],0,"xu"));
  const DoubleVector alpha0  = getDoubleVector(prhs[4]);
  const DoubleVector mu0     = getDoubleVector(prhs[5]);
  const DoubleVector Xr0     = getDoubleVector(prhs[6]);
  const DoubleVector I       = getDoubleVector(prhs[7]);

  // Get the number of samples (n), the number of variables (p), adn
  // the number of coordinate ascent updates (m).
  const Size n = X.nr;
  const Size p = X.nc;
  const Size m = I.n;

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
    Index k = (Index) I.elems[j];

    // Copy the kth column of matrix X.
    copyColumn(X.elems,x,k,n);

    // Perform the update.
    varbvsbinupdate(x,xy.elems[k],xu.elems[k],d.elems[k],u.elems,sa,
		    logodds.elems[k],alpha.elems+k,mu.elems+k,Xr.elems,n);
  }

  // Free the dynamically allocated memory.
  free(x);
}
