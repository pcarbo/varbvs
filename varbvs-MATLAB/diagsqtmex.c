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
// For a description of this C code, see diagsqt.m.
#include "types.h"
#include "doublevectormex.h"
#include "singlematrixmex.h"
#include "misc.h"
#include "diagsq.h"
#include "mex.h"
#include "matrix.h"

// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------
// MEX-file gateway routine.
void mexFunction (int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[]) {

  // Get the input matrix X.
  const SingleMatrix X = getSingleMatrix(prhs[0]);

  // Get the input vector a.
  const DoubleVector a = getDoubleVector(prhs[1]);
    
  // Initialize the output.
  DoubleVector y = createMatlabVector(X.nr,&plhs[0]);
  
  // Compute y = X.^2*a.
  diagsqt(X.elems,a.elems,y.elems,X.nr,X.nc);    
}
