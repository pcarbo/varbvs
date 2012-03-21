#include "varbvsbin.h"
#include "sigmoid.h"
#include "vectorops.h"
#include <math.h>

// Function definitions.
// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in logistic
// regression.
void varbvsbinupdate (const double* x, double xy, double xu, double d, 
		      const double* u, double sa, double logodds, 
		      double* alpha, double* mu, double* Xr, Size n) {

  // Compute the variational estimate of the posterior variance.
  double s = sa/(sa*d + 1);      
  
  // Update the variational estimate of the posterior mean.
  double r = (*alpha) * (*mu);
  *mu = s*(xy + d*r + xu*dot(u,Xr,n)/sum(u,n) - dotscaled(x,Xr,u,n));
    
  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR = (*mu) * (*mu) / s;
  *alpha = sigmoid(logodds + (log(s/sa) + SSR)/2);
  
  // Update Xr = X*r.
  double rnew = (*alpha) * (*mu);
  add(Xr,rnew - r,x,n);
}
