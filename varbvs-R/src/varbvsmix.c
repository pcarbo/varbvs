#include "varbvsmix.h"
#include "sigmoid.h"
#include "vectorops.h"
#include <math.h>

// Function definitions.
// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for the Bayesian variable selection mixture
// model with linear regression.
void varbvsmixupdate (const double* x, double xy, double d, double sigma, 
		      double sa1, double sa2, double logodds, double* alpha, 
		      double* mu1, double* mu2, double* Xr, Size n) {

  // Compute the variational estimates of the posterior variances.
  double s1 = sa1*sigma/(sa1*d + 1);
  double s2 = sa2*sigma/(sa2*d + 1);

  // Update the variational estimates of the posterior means.
  double r = (*alpha) * (*mu1) + (1 - *alpha) * (*mu2);
  *mu1 = s1/sigma * (xy + d*r - dot(x,Xr,n));
  *mu2 = *mu1 * s2/s1;
  
  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR1 = (*mu1) * (*mu1) / s1;
  double SSR2 = (*mu2) * (*mu2) / s2;
  *alpha = sigmoid(logodds + (log(s1*sa2/(s2*sa1)) + SSR1 - SSR2)/2);
  
  // Update Xr = X*r.
  double rnew = (*alpha) * (*mu1) + (1 - *alpha) * (*mu2);
  add(Xr,rnew - r,x,n);
}
