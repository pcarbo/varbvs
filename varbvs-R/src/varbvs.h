#ifndef INCLUDE_VARBVS
#define INCLUDE_VARBVS

#include "types.h"

// Function declarations.
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
void varbvsupdate (const double* x, double xy, double d, double sigma, 
		   double sa, double logodds, double* alpha, double* mu, 
		   double* Xr, Size n);

#endif
