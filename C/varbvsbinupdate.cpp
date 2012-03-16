// For a description of this C++ code, see varbvsbinupdate.m.
#include "common.h"
#include "doublevector.h"
#include "singlematrix.h"
#include <string.h>
#include <math.h>

// These include files have a bunch of definitions to interface C++
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Function declarations.
// -----------------------------------------------------------------
// Execute a single iteration of the coordinate ascent updates.
void varbvsbinupdate (const singlematrix& X, const double* xy, 
		      const double* xu, const double* d, 
		      const doublevector& u, double sa, 
		      const double* logodds, doublevector& alpha, 
		      doublevector& mu, doublevector& Xr, 
		      mwSize m, const double* I);

// Function definitions.
// -----------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {
  const mxArray* ptr = prhs[0];

  // Make sure we have the correct number of inputs and outputs.
  if (nrhs != 8)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 3)
    mexErrMsgTxt("Incorrect number of output arguments");

  // (1.) GET INPUTS.
  // Get the input matrix X.
  if (mxGetClassID(ptr) != mxSINGLE_CLASS)
    mexErrMsgTxt("Input argument X must be SINGLE");
  const singlematrix X = getsinglematrix(ptr);

  // Get the number of samples (n) and number of variables (p).
  const mwSize n = X.nr;
  const mwSize p = X.nc;

  // Get input scalar sa.
  ptr = prhs[1];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != 1)
    mexErrMsgTxt("Input argument SA must be a double precision scalar");
  const double sa = *mxGetPr(ptr);

  // Get input vector logodds.
  ptr = prhs[2];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("LOGODDS must be a double precision vector of length P");
  const doublevector logodds = getdoublevector(ptr);

  // Get the sufficient statistics.
  ptr = prhs[3];
  if (!mxIsStruct(ptr))
    mexErrMsgTxt("Input argument STATS must be a structure array");

  // Get the input vector u.
  const mxArray* field = mxGetField(ptr,0,"u");
  if (!field)
    mexErrMsgTxt("Field U in structure STATS does not exist");
  if (!mxIsDouble(field) || mxGetNumberOfElements(field) != n)
    mexErrMsgTxt("STATS.U must be a double precision vector of length N");
  const doublevector u = getdoublevector(field);

  // Get input vector xy.
  field = mxGetField(ptr,0,"xy");
  if (!field)
    mexErrMsgTxt("Field XY in structure STATS does not exist");
  if (!mxIsDouble(field) || mxGetNumberOfElements(field) != p)
    mexErrMsgTxt("STATS.XY must be a double precision vector of length P");
  const doublevector xy = getdoublevector(field);

  // Get input vector xu.
  field = mxGetField(ptr,0,"xu");
  if (!field)
    mexErrMsgTxt("Field XU in structure STATS does not exist");
  if (!mxIsDouble(field) || mxGetNumberOfElements(field) != p)
    mexErrMsgTxt("STATS.XU must be a double precision vector of length P");
  const doublevector xu = getdoublevector(field);  

  // Get input vector d.
  field = mxGetField(ptr,0,"d");
  if (!field)
    mexErrMsgTxt("Field D in structure STATS does not exist");
  if (!mxIsDouble(field) || mxGetNumberOfElements(field) != p)
    mexErrMsgTxt("STATS.D must be a double precision vector of length P");
  const doublevector d = getdoublevector(field);

  // Get input vector alpha0.
  ptr = prhs[4];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("ALPHA0 must be a double precision vector of length P");  
  const doublevector alpha0 = getdoublevector(ptr);

  // Get input vector mu0.
  ptr = prhs[5];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input MU0 must be a double precision vector of length P");  
  const doublevector mu0 = getdoublevector(ptr);

  // Get the input vector Xr0.
  ptr = prhs[6];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != n)
    mexErrMsgTxt("Input XR0 must be a double precision vector of length N");
  const doublevector Xr0 = getdoublevector(ptr);

  // Get input vector I.
  ptr = prhs[7];
  if (!mxIsDouble(ptr))
    mexErrMsgTxt("Input I must be a double precision vector.");
  const doublevector I = getdoublevector(ptr);

  // Get the number of coordinate ascent updates.
  mwSize m = mxGetNumberOfElements(ptr);

  // (2.) INITIALIZE OUTPUTS.
  doublevector alpha = creatematlabvector(p,plhs[0]);
  doublevector mu    = creatematlabvector(p,plhs[1]);
  doublevector Xr    = creatematlabvector(n,plhs[2]);

  copy(alpha0,alpha);
  copy(mu0,mu);
  copy(Xr0,Xr);

  // (3.) RUN COORDINATE ASCENT UPDATES.
  varbvsbinupdate(X,xy.elems,xu.elems,d.elems,u,sa,logodds.elems,
		  alpha,mu,Xr,m,I.elems);
}

// Execute a single iteration of the coordinate ascent updates.
void varbvsbinupdate (const singlematrix& X, const double* xy, 
		      const double* xu, const double* d, 
		      const doublevector& u, double sa, 
		      const double* logodds, doublevector& alpha, 
		      doublevector& mu, doublevector& Xr, 
		      mwSize m, const double* I) {

  // These variables store some temporary results.
  double  alphak, muk, rk, sk;
  double  SSR;
  mwIndex k;

  // Get the number of samples (n) and the number of variables (p).
  const mwSize n = X.nr;
  const mwSize p = X.nc;

  // This is storage for columns of X.
  doublevector xk = newdoublevector(n);

  // Compute the sum of the entries of vector u.
  double ubar = sum(u);

  // Repeat for each coordinate ascent update.
  for (mwIndex j = 0; j < m; j++) {

    // We need to subtract 1 from the index here because Matlab arrays
    // start at one, but C++ arrays start at zero.
    k = (mwIndex) I[j] - 1;
    if (k < 0 || k >= p)
      mexErrMsgTxt("Input I contains an invalid index");

    // Get the current variational parameters.
    alphak = alpha.elems[k];
    muk    = mu.elems[k]; 
    rk     = alphak * muk;
    
    // Get the kth column of genotype matrix X.
    getcolumn(X,xk,k);
    
    // Compute the variational estimate of the posterior variance.
    sk = sa/(sa*d[k] + 1);      

    // Update the variational estimate of the posterior mean.
    muk = sk * (xy[k] - dot(xk,Xr,u) + xu[k]*dot(u,Xr)/ubar + d[k]*rk);
    
    // Update the variational estimate of the posterior inclusion
    // probability.
    SSR    = muk*muk/sk;
    alphak = sigmoid(logodds[k] + (log(sk/sa) + SSR)/2);
    
    // Update Xr = X*r.
    add(Xr,alphak*muk - rk,xk);

    // Store the updated variational parameters.
    alpha.elems[k] = alphak;
    mu.elems[k]    = muk;
  }

  // Free dynamically allocated memory.
  free(xk);
}
