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
#include "singlematrixmex.h"

// FUNCTION DEFINITIONS
// -----------------------------------------------------------------
// Get a single precision matrix from a MATLAB array. 
SingleMatrix getSingleMatrix (const mxArray* ptr) {
  SingleMatrix result;
  result.nr = mxGetM(ptr);
  result.nc = mxGetN(ptr);
  result.elems = (MatrixElem*) mxGetPr(ptr);
  return result;
}
