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
#ifndef INCLUDE_DOUBLEVECTORMEX
#define INCLUDE_DOUBLEVECTORMEX

#include "types.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// TYPE DEFINITIONS
// -----------------------------------------------------------------
// A vector with double precision entries.
typedef struct {
  Size    n;      // Size of vector.
  double* elems;  // Vector entries.
} DoubleVector;

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Get a double precision vector from a MATLAB array.
DoubleVector getDoubleVector (const mxArray* ptr);

// Create a column vector in a MATLAB array.
DoubleVector createMatlabVector (Size n, mxArray** ptr);

// Copy entries of one vector to another vector.
void copyDoubleVector (const DoubleVector source, DoubleVector dest);

#endif
