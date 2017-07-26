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
#include "diagsq.h"
#include "misc.h"

// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------
// Compute (X.^2)'*a and store the result in vector y.
void diagsq (const MatrixElem* X, const double* a, double* y, 
	     Size m, Size n) {
  double t;

  // Repeat for each column of X.
  for (Index j = 0; j < n; j++, y++) {
    
    // Initialize the jth entry of the result.
    *y = 0;

    // Repeat for each row of X.
    const double* ai = a;
    for (Index i = 0; i < m; i++, X++, ai++) {

      // Add X(i,j)^2 * a(i) to the jth entry of y.
      t   = (double) *X;
      *y += t * t * (*ai);
    }
  }
}

// ---------------------------------------------------------------------
// Compute X.^2*a and store the result in vector y.
void diagsqt (const MatrixElem* X, const double* a, double* y, Size m,
	      Size n) {
  double t;  // An intermediate result.

  // Zero the entries of the result vector.
  setVector(y,m,0);

  // Repeat for each column of X.
  for (Index j = 0; j < n; j++, a++) {

    // Repeat for each row of X.
    double* yi = y;
    for (Index i = 0; i < m; i++, X++, yi++) {
      
      // Add X(i,j)^2 * a(j) to the ith entry of y.
      t    = (double) *X;
      *yi += t * t * (*a);
    }
  }
}
