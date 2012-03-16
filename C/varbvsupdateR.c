// For a description of this C code, see the comments for function
// varbvsupdate in file varbvs.R.
#include "vectorops.h"
#include "varbvs.h"

// These include files have a bunch of definitions to interface C
// routines to R.
#include <R.h>
#include <Rinternals.h>

void varbvsupdateR (const int* np, const int* mp, const const double* X, 
		    const double* sigmap, const double* sap, 
		    const double* logodds, const double* xy, 
		    const double* d, const double* alpha0, 
		    const double* my0, const double* Xr0, const int* I, 
		    double* alpha, double* mu, double* Xr, double* x) {

  // Get the number of samples (n) and the number of coordinate ascent
  // updates (m).
  const int n = *np;
  const int m = *mp;

  // Get input scalars sigma and sa.
  const double sigma = *sigmap;
  const double sa    = *sap;

  // Run the coordinate ascent updates.
  for (int j = 0; j < m; j++) {
    int k = I[j];

    // Get the kth column of matrix X.
    getColumn(X,x,k,n);

    // Perform the update.
    varbvsupdate(x,xy[k],d[k],sigma,sa,logodds[k],alpha+k,mu+k,Xr,n);
  }
}
