// For a description of this C code, see the help for function
// 'varbvsbinupdate' in R.
#include "vectorops.h"
#include "varbvsbin.h"

// These include files have a bunch of definitions to interface C
// routines to R.
#include <R.h>
#include <Rinternals.h>

void varbvsupdateR (const int* np, const int* mp, const double* X, 
		    const double* sigmap, const double* sap, 
		    const double* logodds, const double* xy, 
		    const double* xu,
		    const double* d, double* alpha, double* mu, 
		    double* Xr, const int* S) {

  // // Get the number of samples (n) and the number of coordinate ascent
  // // updates (m).
  // const int n = *np;
  // const int m = *mp;

  // // Get input scalars sigma and sa.
  // const double sigma = *sigmap;
  // const double sa    = *sap;

  // Run the coordinate ascent updates.
  for (int j = 0; j < m; j++) {
    int k = S[j];

    // Get the kth column of matrix X.
    const double* xk = getColumn(X,k,n);

    // Perform the update.
    varbvsbinupdate(xk,xy[k],xu[k],d[k], 
		     const double* u, double sa, double logodds, 
		      double* alpha, double* mu, double* Xr, Size n);
  //   varbvsupdate(xk,xy[k],d[k],sigma,sa,logodds[k],alpha+k,mu+k,Xr,n);
  }
}
