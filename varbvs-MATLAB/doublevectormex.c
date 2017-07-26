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
#include "doublevectormex.h"
#include "misc.h"

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// Get a double precision vector from a MATLAB array.
DoubleVector getDoubleVector (const mxArray* ptr) {
  DoubleVector result;
  result.n = mxGetNumberOfElements(ptr);
  result.elems = mxGetPr(ptr);
  return result;
}

// Create a column vector in MATLAB array.
DoubleVector createMatlabVector (Size n, mxArray** ptr) {
  *ptr = mxCreateDoubleMatrix(n,1,mxREAL);
  return getDoubleVector(*ptr);
}

// Copy entries of one vector to another vector.
void copyDoubleVector (const DoubleVector source, DoubleVector dest) {
  copy(source.elems,dest.elems,dest.n);
}
