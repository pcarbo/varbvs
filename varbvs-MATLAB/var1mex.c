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
// For a description of this C code, see var1.m.
#include "types.h"
#include "misc.h"
#include "doublevectormex.h"
#include "singlematrixmex.h"
#include "mex.h"
#include "matrix.h"

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Compute the sample variance.
double computeVariance (const double* x, Size n);

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // Get the input matrix X.
  const SingleMatrix X = getSingleMatrix(prhs[0]);

  // Get the number of samples (n) and the number of variables (p).
  const Size n = X.nr;
  const Size p = X.nc;

  // Initialize the output.
  DoubleVector y = createMatlabVector(p,&plhs[0]);

  // This is storage for a column of matrix X.
  double* x = malloc(sizeof(double)*n);

  // Repeat for each column of X.
  for (Index k = 0; k < p; k++) {

    // Get the kth column of X.
    copyColumn(X.elems,x,k,n);

    // Compute the sample variance.
    y.elems[k] = computeVariance(x,n);
  }

  // Free the dynamically allocated memory.
  free(x);
}

// -----------------------------------------------------------------
// Compute the sample variance.
double computeVariance (const double* x, Size n) {
  double y = 0;  // The return value.
  double t;      // Intermediate result.

  // Compute the sample mean.
  double mu = sum(x,n) / n;

  // Repeat for each entry of x.
  for (Index i = 0; i < n; i++) {
    t  = x[i] - mu;
    y += t*t;
  }

  // Divide by the number of entries.
  y /= n;

  return y;
}
