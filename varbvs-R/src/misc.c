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
#include "misc.h"
#include <string.h>
#include <math.h>

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// Computes log(1 + exp(x)).
double logpexp (double x) {
  return (x >= 16) * x + (x < 16)  * log(1 + exp(x));
}

// Return the sigmoid function at x.
double sigmoid (double x) {
  return 1/(1 + exp(-x));
}

// Return the logarithm of the sigmoid function at x.
double logsigmoid (double x) {
  return -logpexp(-x);
}

// This function takes as input an array of unnormalized
// log-probabilities "logw" and returns normalized probabilities "w"
// such that the sum is equal to 1.
void normalizelogweights (const double* logw, double* w, Size n) {

  // Guard against underflow or overflow by adjusting the
  // log-probabilities so that the largest probability is 1.
  double c = max(logw,n);
  for (Index i = 0; i < n; i++)
    w[i] = exp(logw[i] - c);

  // Normalize the probabilities.
  double r = sum(w,n);
  for (Index i = 0; i < n; i++)
    w[i] /= r;
}
  
// Copy entries of one vector to another vector.
void copy (const double* source, double* dest, Size n) {
  memcpy(dest,source,sizeof(double)*n);
}

// Get a pointer to column j of matrix X.
const MatrixElem* getConstColumn (const MatrixElem* X, Index j, Size n) {
  return X + n*j;
}

// Get a pointer to column j of matrix X.
MatrixElem* getColumn (MatrixElem* X, Index j, Size n) {
  return X + n*j;
}

// Copy column j of matrix X.
void copyColumn (const MatrixElem* X, double* y, Index j, Size n) {
  const MatrixElem* xij = getConstColumn(X,j,n);
  for (Index i = 0; i < n; i++, xij++, y++)
    *y = (double) *xij;
}

// Set the entries of the vector to a.
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
  double y = x[0];
  for (Index i = 1; i < n; i++)
    y = (x[i] > y) * x[i] + (x[i] <= y) * y;
  return y;
}

// Add a*x to vector y, and store the result in y.
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

// Compute X'*y and return the result in a.
void matrixvec (const double* X, const double* y, double* a, 
		Size nr, Size nc) {
  Index i = 0;
  for (Index c = 0; c < nc; c++) {
    a[c] = 0;
    for (Index r = 0; r < nr; r++, i++)
      a[c] += X[i] * y[r];
  }
}
