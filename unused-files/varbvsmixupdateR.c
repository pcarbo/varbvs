#include "vectorops.h"
#include "varbvsmix.h"

// These include files have a bunch of definitions to interface C
// routines to R.
#include <R.h>
#include <Rinternals.h>

void varbvsmixupdateR (const int* np, const int* mp, const double* X, 
		       const double* sigmap, const double* sa1p,
		       const double* sa2p, const double* logodds, 
		       const double* xy, const double* d, double* alpha, 
		       double* mu1, double* mu2, double* Xr, const int* S) {

  // Get the number of samples (n) and the number of coordinate ascent
  // updates (m).
  const int n = *np;
  const int m = *mp;

  // Get input scalars sigma, sa1 and sa2.
  const double sigma = *sigmap;
  const double sa1   = *sa1p;
  const double sa2   = *sa2p;

  // Run the coordinate ascent updates.
  for (int j = 0; j < m; j++) {
    int k = S[j];

    // Get the kth column of matrix X.
    const double* xk = getColumn(X,k,n);

    // Perform the update.
    varbvsmixupdate(xk,xy[k],d[k],sigma,sa1,sa2,logodds[k],
		    alpha+k,mu1+k,mu2+k,Xr,n);
  }
}
