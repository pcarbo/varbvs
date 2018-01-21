// Part of the varbvs package, https://github.com/pcarbo/varbvs
//
// Copyright (C) 2012-2018, Peter Carbonetto
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
#ifndef INCLUDE_VARBVS
#define INCLUDE_VARBVS

#include "types.h"

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in linear
// regression. The inputs are as follows: x is a single column of the
// data matrix X; xy is the corresponding entry of matrix-vector
// product X'*y; d is the corresponding entry on the diagonal of X'*X;
// sigma, sa and logodds specify the hyperparameters; n is the number
// of samples; alpha and mu are the variational parameters that will
// be updated; and Xr is the matrix-vector product X'*r which will be
// updated to reflect the change to alpha and mu.
void varbvsnormupdate (const double* x, double xy, double d, double sigma, 
		       double sa, double logodds, double* alpha, double* mu, 
		       double* Xr, Size n);

// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in logistic
// regression. The inputs are as follows: x is a single column of the
// data matrix X; d is a vector defined in the Bayesian Analysis paper
// that is used to specify the variational approximation to the
// logistic regression factors; xy is the corresponding entry of
// matrix-vector product X'*yhat; xd is the corresponding entry of
// matrix-vector product X'*d; xdx is the corresponding entry on the
// diagonal of X'*D*X; sa and logodds specify the hyperparameters; n is
// the number of samples; alpha and mu are the variational parameters
// that will be updated; and Xr is the matrix-vector product X*r
// which will be updated to reflect the changes to alpha and mu.
void varbvsbinupdate (const double* x, double xy, double xd, double xdx, 
		      const double* d, double sa, double logodds, 
		      double* alpha, double* mu, double* Xr, Size n);

// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in logistic
// regression, allowing for covariates. The inputs are as follows: x
// is a single column of the data matrix X; d is a vector defined in
// the Bayesian Analysis paper that is used to specify the variational
// approximation to the logistic regression factors; xy is the
// corresponding entry of the matrix-vector product X'*yhat; xdx is
// the corresponding entry on the diagonal of X'*D*X; sa and logodds
// specify the hyperparameters; n is the number of samples; m is the
// number of covariates; dxr is the n x m matrix D*Z*R', where R is
// the upper triangular Choleksy factor of S; alpha and mu are the
// variational parameters that will be updated; and Xr is the
// matrix-vector product X*r which will be udpated to reflect the
// changes to alpha and mu.
void varbvsbinzupdate (const double* x, double xy, double xdx, 
		       const double* d, const double* dzr, double sa, 
		       double logodds, double* alpha, double* mu, 
		       double* Xr, double* a, double* b, Size n, Size m);

// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for the linear regression model with
// mixture-of-normal priors. The inputs are as follows: x is a single
// column of the data matrix X; xy is the corresponding entry of
// matrix-vector product X'*y; d is the corresponding entry on the
// diagonal of X'*X; sigma, sa and q specify the hyperparameters (the
// residual variance, the vector of mixture variances, and the vector
// mixture weights, respectively); n is the number of samples; k is
// the number of mixture components (including the "spike"); alpha and
// mu are the variational parameters that will be updated; Xr is the
// matrix-vector product X'*r which will be updated to reflect the
// change to alpha and mu; inputs s and logw are vectors of
// length k that will be used to store intermediate calculations; and
// eps is the floating-point precision value (e.g., type "help eps" in
// MATLAB).
//
// Note that the variational parameters alpha and mu are stored as k x
// p matrices, where p is the number of variables. The variational
// parameters are stored in this way so that that the co-ordinate
// ascent updates are easier to implement.
void varbvsmixupdate (const double* x, double xy, double d, double sigma, 
		      const double* sa, const double* q, double* alpha,
		      double* mu, double* Xr, double* s, double* logw,
		      Size n, Size k, double eps);

#endif
