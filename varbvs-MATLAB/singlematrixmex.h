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
#ifndef INCLUDE_SINGLEMATRIXMEX
#define INCLUDE_SINGLEMATRIXMEX

#include "types.h"
#include "mex.h"
#include "matrix.h"

// TYPE DEFINIITIONS
// -----------------------------------------------------------------
// A dense matrix with single precision floating point entries.
typedef struct {
  Size        nr;     // Number of rows.
  Size        nc;     // Number of columns.
  MatrixElem* elems;  // Entries of matrix.
} SingleMatrix;

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Get a single precision matrix from a MATLAB array. 
SingleMatrix getSingleMatrix (const mxArray* ptr);

#endif
