#ifndef INCLUDE_VARBVS
#define INCLUDE_VARBVS

#include "types.h"

// Function declarations.
// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for Bayesian variable selection in logistic
// regression. The inputs are as follows: x is a single column of the
// data matrix X; u is a vector defined in the Bayesian Analysis paper
// that is used to decribe the variational approximation to the
// logistic regression factors; xy is the corresponding entry of
// matrix-vector product X'*y; xu is the corresponding entry of
// matrix-vector product X'*u; d is the corresponding entry on the
// diagonal of X'*X; sa and logodds specify the hyperparameters; n is
// the number of samples; alpha and mu are the variational parameters
// that will be updated; and Xr is the matrix-vector product X'*r
// which will be updated to reflect the change to alpha and mu.
void varbvsbinupdate (const double* x, double xy, double xu, double d, 
		      const double* u, double sa, double logodds, 
		      double* alpha, double* mu, double* Xr, Size n);

#endif
