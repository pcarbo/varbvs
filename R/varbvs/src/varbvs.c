#include "varbvs.h"
#include "sigmoid.h"
#include "vectorops.h"
#include <math.h>

// Function definitions.
// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in linear
// regression.
void varbvsupdate (const double* x, double xy, double d, double sigma, 
		   double sa, double logodds, double* alpha, double* mu, 
		   double* Xr, Size n) {

  // Compute the variational estimate of the posterior variance.
  double s = sa*sigma/(sa*d + 1);
  
  // Update the variational estimate of the posterior mean.
  double r = (*alpha) * (*mu);
  *mu = s/sigma * (xy + d*r - dot(x,Xr,n));
  
  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR = (*mu) * (*mu) / s;
  *alpha = sigmoid(logodds + (log(s/(sa*sigma)) + SSR)/2);
  
  // Update Xr = X*r.
  double rnew = (*alpha) * (*mu);
  add(Xr,rnew - r,x,n);
}
