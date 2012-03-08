// For a description of this C++ code, see varbvsupdate.m.
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
void varbvsupdate (const singlematrix& X, const double* xy, 
		   const double* d, double sigma, double sb, 
		   const double* logodds, doublevector& alpha, 
		   doublevector& mu, doublevector& Xr, mwSize m, 
		   const double* snps);

// -----------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {
  const mxArray* ptr = prhs[0];

  // Make sure we have the correct number of inputs and outputs.
  if (nrhs != 10)
    mexErrMsgTxt("Incorrect number of input arguments");
  if (nlhs != 3)
    mexErrMsgTxt("Incorrect number of output arguments");

  // (1.) GET INPUTS.
  // Get the input matrix X.
  if (mxGetClassID(ptr) != mxSINGLE_CLASS)
    mexErrMsgTxt("Input argument X must be SINGLE");
  const singlematrix X = getsinglematrix(ptr);

  // Get the size of the data set (n) and the number of SNPs (p).
  const int n = X.nr;
  const int p = X.nc;

  // Get input scalar sigma.
  ptr = prhs[1];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != 1)
    mexErrMsgTxt("Input argument SIGMA must be a double precision scalar");
  const double sigma = *mxGetPr(ptr);

  // Get input scalar sb.
  ptr = prhs[2];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != 1)
    mexErrMsgTxt("Input argument SB must be a double precision scalar");
  const double sb = *mxGetPr(ptr);

  // Get input vector logodds.
  ptr = prhs[3];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument LOGODDS must be a double precision \
vector of length P");
  const doublevector logodds = getdoublevector(ptr);

  // Get input vector xy.
  ptr = prhs[4];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument XY must be a double precision vector \
of length P");
  const doublevector xy = getdoublevector(ptr);

  // Get input vector d.
  ptr = prhs[5];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument D must be a double precision vector \
of length P");
  const doublevector d = getdoublevector(ptr);

  // Get input vector alpha0.
  ptr = prhs[6];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument A0 must be a double precision \
vector of length P");  
  const doublevector alpha0 = getdoublevector(ptr);

  // Get input vector mu0.
  ptr = prhs[7];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument MU0 must be a double precision vector \
of length P");  
  const doublevector mu0 = getdoublevector(ptr);

  // Get the input vector Xr0.
  ptr = prhs[8];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != n)
    mexErrMsgTxt("Input argument XR0 must be a double precision vector \
of length N");
  const doublevector Xr0 = getdoublevector(ptr);

  // Get input vector snps.
  ptr = prhs[9];
  if (!mxIsDouble(ptr))
    mexErrMsgTxt("Input argument KS must be a double precision vector.");
  const doublevector snps = getdoublevector(ptr);

  // Get the number of coordinate ascent updates.
  mwSize m = mxGetNumberOfElements(ptr);

  // (2.) INITIALIZE OUTPUTS.
  doublevector alpha = creatematlabvector(p,plhs[0]);
  doublevector mu    = creatematlabvector(p,plhs[1]);
  doublevector Xr    = creatematlabvector(n,plhs[2]);

  copy(alpha0,alpha);
  copy(mu0,mu);
  copy(Xr0,Xr);

  // (4.) RUN COORDINATE ASCENT UPDATES.
  multisnpupdate(X,xy.elems,d.elems,sigma,sb,logodds.elems,alpha,mu,
		 Xr,m,snps.elems);
}

// Execute a single iteration of the coordinate ascent updates.
void multisnpupdate (const singlematrix& X, const double* xy, 
		     const double* d, double sigma, double sb, 
		     const double* logodds, doublevector& alpha, 
		     doublevector& mu, doublevector& Xr, mwSize m, 
		     const double* snps) {

  // These variables store some temporary results.
  double  alphak, muk, rk, sk;
  double  SSR;
  mwIndex k;

  // Get the size of the data set (n) and the number of SNPs (p).
  const mwSize n = X.nr;
  const mwSize p = X.nc;

  // This is storage for columns of X.
  doublevector xk = newdoublevector(n);

  // Repeat for each coordinate ascent update.
  for (mwIndex j = 0; j < m; j++) {

    // We need to subtract 1 from the SNP index here because Matlab
    // arrays start at one, but C++ arrays start at zero.
    k = (mwIndex) snps[j] - 1;

    // Get the current variational parameters.
    alphak = alpha.elems[k];
    muk    = mu.elems[k]; 
    rk     = alphak * muk;

    // Get the kth column of genotype matrix X.
    getcolumn(X,xk,k);

    // Compute variational parameter S.
    sk = sb*sigma/(sb*d[k] + 1);

    // Update variational parameter MU.
    muk = sk/sigma * (xy[k] + d[k]*rk - dot(xk,Xr));

    // Update variational parameter ALPHA.
    SSR    = muk*muk/sk;
    alphak = sigmoid(logodds[k] + (log(sk/(sb*sigma)) + SSR)/2);
    
    // Update Xr = X*r.
    add(Xr,alphak*muk - rk,xk);

    // Store the updated parameters.
    alpha.elems[k] = alphak;
    mu.elems[k]    = muk;
  }

  // Free dynamically allocated memory.
  free(xk);
}
