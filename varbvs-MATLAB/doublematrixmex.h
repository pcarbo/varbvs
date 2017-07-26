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
#ifndef INCLUDE_DOUBLEMATRIXMEX
#define INCLUDE_DOUBLEMATRIXMEX

#include "types.h"
#include "mex.h"
#include "matrix.h"

// TYPE DEFINITIONS
// -----------------------------------------------------------------
// A dense matrix with single precision entries.
typedef struct {
  Size    nr;     // Number of rows.
  Size    nc;     // Number of columns.
  double* elems;  // Entries of matrix.
} DoubleMatrix;

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Get a double precision floating point matrix from a MATLAB array. 
DoubleMatrix getDoubleMatrix (const mxArray* ptr);

// Create an m x n matrix in a MATLAB array.
DoubleMatrix createMatlabDoubleMatrix (Size m, Size n, mxArray** ptr);

// Copy all the entries of one matrix to another matrix.
void copyDoubleMatrix (const DoubleMatrix source, DoubleMatrix dest);

// Get a pointer to column j of matrix X. The entries in a single
// column of the matrix are assumed to be stored consecutively in
// memory. Input n is the number of rows in the matrix.
const double* getConstDoubleColumn (const double* X, Index j, Size n);
  
// Get a pointer to column j of matrix X. The entries in a single
// column of the matrix are assumed to be stored consecutively in
// memory. Input n is the number of rows in the matrix.
double* getDoubleColumn (double* X, Index j, Size n);

// Copy column j of matrix X. The entries in a single column of the
// matrix are assumed to be stored consecutively in memory. Input n
// is the number of rows in the matrix.
void copyDoubleColumn (const double* X, double* y, Index j, Size n);

#endif
