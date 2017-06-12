// Part of the varbvs package, https://github.com/pcarbo/varbvs
//
// Copyright (C) 2012-2017, Peter Carbonetto
//
// This program is free software: you can redistribute it under the
// terms of the GNU General Public License; either version 3 of the
// License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANY; without even the implied warranty of
// MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
// General Public License for more details.
//
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
  double s = sigma*sa/(sa*d + 1);
  
  // Update the variational estimate of the posterior mean.
  double r = (*alpha) * (*mu);
  *mu = s/sigma * (xy + d*r - dot(x,Xr,n));
  
  // Update the variational estimate of the posterior inclusion
  // probability.
  double SSR = (*mu) * (*mu) / s;
  *alpha = sigmoid(logodds + (log(s/(sigma*sa)) + SSR)/2);
  
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

// ---------------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the variational 
// lower bound for the linear regression model with mixture-of-normal
// priors.
void varbvsmixupdate (const double* x, double xy, double d, double sigma, 
		      const double* sa, const double* q, double* alpha,
		      double* mu, double* Xr, double* s, double* logw,
		      Size n, Size k, double eps) {

  // The mean and variance corresponding to the first mixture
  // component, the "spike", should always be zero.
  mu[0]   = 0;
  s[0]    = 0;
  logw[0] = log(q[0] + eps);
  
  // Compute the variance of the regression coefficient conditioned on
  // being drawn from each of the mixture components.
  for (Index i = 1; i < k; i++)
    s[i] = sigma*sa[i]/(sa[i]*d + 1);
  
  // Update the variational estimate of the posterior mean for each
  // mixture component.
  double r = dot(alpha,mu,k);
  double t = xy + d*r - dot(x,Xr,n);
  for (Index i = 1; i < k; i++)
    mu[i] = s[i]/sigma*t; 

  // Update the assignment probabilities for all of the mixture
  // components.
  double SSR;
  for (Index i = 1; i < k; i++) {
    SSR     = mu[i]*mu[i]/s[i];
    logw[i] = log(q[i] + eps) + (log(s[i]/(sigma*sa[i])) + SSR)/2;
  }
  normalizelogweights(logw,alpha,k);
  
  // Update Xr = X*r.
  double rnew = dot(alpha,mu,k);
  add(Xr,rnew - r,x,n);
}
