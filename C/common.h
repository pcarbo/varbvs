#ifndef INCLUDE_COMMON
#define INCLUDE_COMMON

// This include file has a bunch of definitions to interface C++
// routines to MATLAB.
#include "matrix.h"

// Function declarations.
// -----------------------------------------------------------------
// Computes log(1 + exp(x)) in a numerically stable manner.
double logpexp (double x);

// Returns the sigmoid function at x.
double sigmoid (double x);

// Returns the logarithm of the sigmoid function at x. Computation is
// performed in a numerically stable manner.
double logsigmoid (double x);

// Return the sum of the n elements in array x.
double sum (const double* x, mwSize n);

// Return the maximum of the n elements in array x.
double max (const double* x, mwSize n);

#endif
