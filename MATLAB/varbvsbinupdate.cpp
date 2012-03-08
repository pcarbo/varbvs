// For a description of this C++ code, see multisnpbinupdate.m.
#include "mex.h"
#include "matrix.h"
#include "common.h"
#include "doublevector.h"
#include "singlematrix.h"
#include <string.h>
#include <math.h>

// Function declarations.
// -----------------------------------------------------------------
// Execute a single iteration of the coordinate ascent updates.
void multisnpbinupdate (const singlematrix& X, const double* xy, 
			const double* xu, const double* d, 
			const doublevector& u, double sb, 
			const double* logodds, doublevector& alpha, 
			doublevector& mu, doublevector& Xr, 
			mwSize m, const double* snps);

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

  // Get the size of the data set (n) and the number of SNPs (p).
  const mwSize n = X.nr;
  const mwSize p = X.nc;

  // Get input scalar sb.
  ptr = prhs[1];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != 1)
    mexErrMsgTxt("Input argument SB must be a double precision scalar");
  const double sb = *mxGetPr(ptr);

  // Get input vector logodds.
  ptr = prhs[2];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument LOGODDS must be a double precision \
vector of length P");
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
    mexErrMsgTxt("Input argument STATS.U must be a double precision \
vector of length N");
  const doublevector u = getdoublevector(field);

  // Get input vector xy.
  field = mxGetField(ptr,0,"xy");
  if (!mxIsDouble(field) || mxGetNumberOfElements(field) != p)
    mexErrMsgTxt("Input argument STATS.XY must be a double precision vector \
of length P");
  const doublevector xy = getdoublevector(field);

  // Get input vector xu.
  field = mxGetField(ptr,0,"xu");
  if (!mxIsDouble(field) || mxGetNumberOfElements(field) != p)
    mexErrMsgTxt("Input argument STATS.XU must be a double precision vector \
of length P");
  const doublevector xu = getdoublevector(field);  

  // Get input vector d.
  field = mxGetField(ptr,0,"d");
  if (!mxIsDouble(field) || mxGetNumberOfElements(field) != p)
    mexErrMsgTxt("Input argument STATS.D must be a double precision vector \
of length P");
  const doublevector d = getdoublevector(field);

  // Get input vector alpha0.
  ptr = prhs[4];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument ALPHA0 must be a double precision \
vector of length P");  
  const doublevector alpha0 = getdoublevector(ptr);

  // Get input vector mu0.
  ptr = prhs[5];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != p)
    mexErrMsgTxt("Input argument MU0 must be a double precision vector \
of length P");  
  const doublevector mu0 = getdoublevector(ptr);

  // Get the input vector Xr0.
  ptr = prhs[6];
  if (!mxIsDouble(ptr) || mxGetNumberOfElements(ptr) != n)
    mexErrMsgTxt("Input argument XR0 must be a double precision vector \
of length N");
  const doublevector Xr0 = getdoublevector(ptr);

  // Get input vector snps.
  ptr = prhs[7];
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

  // (3.) RUN COORDINATE ASCENT UPDATES.
  multisnpbinupdate(X,xy.elems,xu.elems,d.elems,u,sb,logodds.elems,
		    alpha,mu,Xr,m,snps.elems);
}

// Execute a single iteration of the coordinate ascent updates.
void multisnpbinupdate (const singlematrix& X, const double* xy, 
			const double* xu, const double* d, 
			const doublevector& u, double sb, 
			const double* logodds, doublevector& alpha, 
			doublevector& mu, doublevector& Xr, 
			mwSize m, const double* snps) {

  // These variables store some temporary results.
  double  alphak, muk, rk, sk;
  double  SSR;
  mwIndex k;

  // Get the size of the data set (n) and the number of SNPs (p).
  const mwSize n = X.nr;
  const mwSize p = X.nc;

  // This is storage for columns of X.
  doublevector xk = newdoublevector(n);

  // Compute the sum of the entries of vector u.
  double ubar = sum(u);

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
    sk = sb/(sb*d[k] + 1);      

    // Update variational parameter MU.
    muk = sk * (xy[k] - dot(xk,Xr,u) + xu[k]*dot(u,Xr)/ubar + d[k]*rk);
    
    // Update variational parameter ALPHA.
    SSR    = muk*muk/sk;
    alphak = sigmoid(logodds[k] + (log(sk/sb) + SSR)/2);
    
    // Update Xr = X*r.
    add(Xr,alphak*muk - rk,xk);

    // Store the updated parameters.
    alpha.elems[k] = alphak;
    mu.elems[k]    = muk;
  }

  // Free dynamically allocated memory.
  free(xk);
}
