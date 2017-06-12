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
#include "doublematrixmex.h"
#include "misc.h"

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// Get a double precision matrix from a MATLAB array. 
DoubleMatrix getDoubleMatrix (const mxArray* ptr) {
  DoubleMatrix result;
  result.nr    = mxGetM(ptr);
  result.nc    = mxGetN(ptr);
  result.elems = (double*) mxGetPr(ptr);
  return result;
}

// Create an m x n matrix in a MATLAB array.
DoubleMatrix createMatlabDoubleMatrix (Size m, Size n, mxArray** ptr) {
  *ptr = mxCreateDoubleMatrix(m,n,mxREAL);
  return getDoubleMatrix(*ptr);
}

// Copy all the entries of one matrix to another matrix.
void copyDoubleMatrix (const DoubleMatrix source, DoubleMatrix dest) {
  copy(source.elems,dest.elems,dest.nr*dest.nc);
}

// Get a pointer to column j of matrix X.
const double* getConstDoubleColumn (const double* X, Index j, Size n) {
  return X + n*j;
}

// Get a pointer to column j of n x m matrix X.
double* getDoubleColumn (double* X, Index j, Size n) {
  return X + n*j;
}

// Copy column j of matrix X.
void copyDoubleColumn (const double* X, double* y, Index j, Size n) {
  const double* xij = getConstDoubleColumn(X,j,n);
  for (Index i = 0; i < n; i++, xij++, y++)
    *y = (double) *xij;
}

