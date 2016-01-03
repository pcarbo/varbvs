#ifndef INCLUDE_VARBVSMIX
#define INCLUDE_VARBVSMIX

#include "types.h"

// Function declarations.
// -----------------------------------------------------------------
// Execute a single coordinate ascent update to maximize the
// variational lower bound for the Bayesian variable selection mixture
// model with linear regression.
void varbvsmixupdate (const double* x, double xy, double d, double sigma, 
		      double sa1, double sa2, double logodds, double* alpha, 
		      double* mu1, double* mu2, double* Xr, Size n);

#endif
