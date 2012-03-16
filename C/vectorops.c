#include "vectorops.h"
#include <string.h>

// Function definitions.
// -----------------------------------------------------------------
// Copy entries of one vector to another vector.
void copy (const double* source, double* dest, Size n) {
  memcpy(dest,source,sizeof(double)*n);
}

// Copy column j of matrix X.
void getColumn (const MatrixElem* X, double* y, Index j, Size n) {
  const MatrixElem* xij = X + n*j;  // Pointer to (i,j) entry of X.
  for (Index i = 0; i < n; i++, xij++, y++)
    *y = (double) *xij;
}

// Set the entries of the vector to the same value.
void setVector (double* x, Size n, double a) {
  for (Index i = 0; i < n; i++)
    x[i] = a;
}

// Return the sum of the entries of the vector.
double sum (const double* x, Size n) {
  double y = 0;
  for (Index i = 0; i < n; i++)
    y += x[i];
  return y;
}

// Return the largest entry in the vector.
double max (const double* x, Size n) {
  double y = x[0]; // The return value.
  double xi;
  for (Index i = 1; i < n; i++) {
    xi = x[i];
    y  = (xi > y) * xi + (xi <= y) * y;
  }
  return y;
}

// Add a to all the entries in vector x, and store the result in vector y.
void add (double* y, double a, const double* x, Size n) {
  for (Index i = 0; i < n; i++)
    y[i] += a * x[i];
}

// Return the dot product of two vectors.
double dot (const double* x, const double* y, Size n) {
  double z = 0;
  for (Index i = 0; i < n; i++)
    z += x[i] * y[i];
  return z;
}

// Compute x'*D*y, the dot product of x and y scaled by D = diag(d).
double dotscaled (const double* x, const double* y, const double* d, 
		  Size n) {
  double z = 0;
  for (Index i = 0; i < n; i++)
    z += x[i] * y[i] * d[i];
  return z;
}
