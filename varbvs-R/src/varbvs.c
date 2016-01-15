#include "varbvs.h"
#include "misc.h"
#include <math.h>

// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the variational 
// lower bound for Bayesian variable selection in linear regression.
void varbvsnormupdate (const double* x, double xy, double d, double sigma, 
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

// ---------------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the variational 
// lower bound for Bayesian variable selection in logistic regression.
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

// ---------------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the variational 
// lower bound for Bayesian variable selection in logistic regression, 
// allowing for covariates.
void varbvsbinzupdate (const double* x, double xy, double xdx, 
		       const double* d, const double* dzr, double sa, 
		       double logodds, double* alpha, double* mu, 
		       double* Xr, double* a, double* b, Size n, Size m) {

  // Compute the variational estimate of the posterior variance.
  double s = sa/(sa*xdx + 1);

  // Update the variational estimate of the posterior mean.
  double r = (*alpha) * (*mu);
  matrixvec(dzr,x,a,n,m);
  matrixvec(dzr,Xr,b,n,m);
  *mu = s*(xy + xdx*r + dot(a,b,m) - dotscaled(x,Xr,d,n));

  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR = (*mu) * (*mu) / s;
  *alpha = sigmoid(logodds + (log(s/sa) + SSR)/2);

  // Update Xr = X*r.
  double rnew = (*alpha) * (*mu);
  add(Xr,rnew - r,x,n);
}
