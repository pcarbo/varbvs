#include "common.h"
#include <math.h>

// Function definitions.
// -----------------------------------------------------------------
// Computes log(1 + exp(x)) in a numerically stable manner.
double logpexp (double x) {
  return (x >= 8) * x + (x < 8) * log(1 + exp(x));
}

// Returns the sigmoid function at x.
double sigmoid (double x) {
  return 1/(1 + exp(-x));
}

// Returns the logarithm of the sigmoid function at x. Computation is
// performed in a numerically stable manner.
double logsigmoid (double x) {
  return -logpexp(-x);
}

// Return the sum of the n elements in array x.
double sum (const double* x, mwSize n) {
  double y = 0;  // The return value.

  // Repeat for each element in the array.
  for (mwIndex i = 0; i < n; i++)
    y += x[i];

  return y;
}

// Return the maximum of the n elements in array x.
double max (const double* x, mwSize n) {
  double y = x[0];  // The return value.
  double xi;

  // Repeat for each elements of the array.
  for (mwIndex i = 1; i < n; i++) {
    xi = x[i];
    y  = (xi > y) * xi + (xi <= y) * y;
  }

  return y;
}
