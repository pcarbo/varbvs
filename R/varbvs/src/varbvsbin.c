#include "varbvsbin.h"
#include "sigmoid.h"
#include "vectorops.h"
#include <math.h>

// Function definitions.
// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in logistic
// regression.
void varbvsbinupdate (const double* x, double xy, double xd, double xdx, 
		      const double* d, double sa, double logodds, 
		      double* alpha, double* mu, double* Xr, Size n) {

  // Compute the variational estimate of the posterior variance.
  double s = sa/(sa*xdx + 1);      
  
  // Update the variational estimate of the posterior mean.
  double r = (*alpha) * (*mu);
  *mu = s*(xy + xdx*r + xd*dot(d,Xr,n)/sum(d,n) - dotscaled(x,Xr,d,n));
    
  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR = (*mu) * (*mu) / s;
  *alpha = sigmoid(logodds + (log(s/sa) + SSR)/2);
  
  // Update Xr = X*r.
  double rnew = (*alpha) * (*mu);
  add(Xr,rnew - r,x,n);
}
