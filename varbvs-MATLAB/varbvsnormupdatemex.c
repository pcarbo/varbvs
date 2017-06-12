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
// For a description of this C code, see varbvsnormupdate.m.
#include "types.h"
#include "misc.h"
#include "doublevectormex.h"
#include "singlematrixmex.h"
#include "varbvs.h"
#include "mex.h"
#include "matrix.h"

void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // (1) GET INPUTS
  // --------------
  const SingleMatrix X       = getSingleMatrix(prhs[0]);
  const double       sigma   = *mxGetPr(prhs[1]);
  const double       sa      = *mxGetPr(prhs[2]);
  const DoubleVector logodds = getDoubleVector(prhs[3]);
  const DoubleVector xy      = getDoubleVector(prhs[4]);
  const DoubleVector d       = getDoubleVector(prhs[5]);
  const DoubleVector alpha0  = getDoubleVector(prhs[6]);
  const DoubleVector mu0     = getDoubleVector(prhs[7]);
  const DoubleVector Xr0     = getDoubleVector(prhs[8]);
  const DoubleVector I       = getDoubleVector(prhs[9]);

  // Get the number of samples (n), the number of variables (p), and
  // the number of coordinate ascent updates (m).
  const Size n = X.nr;
  const Size p = X.nc;
  const Size m = I.n;

  // (2) INITIALIZE OUTPUTS
  // ----------------------
  DoubleVector alpha = createMatlabVector(p,&plhs[0]);
  DoubleVector mu    = createMatlabVector(p,&plhs[1]);
  DoubleVector Xr    = createMatlabVector(n,&plhs[2]);

  copyDoubleVector(alpha0,alpha);
  copyDoubleVector(mu0,mu);
  copyDoubleVector(Xr0,Xr);

  // This is storage for a column of matrix X.
  double* x = malloc(sizeof(double)*n);

  // (3) RUN COORDINATE ASCENT UPDATES
  // ---------------------------------
  // Repeat for each coordinate ascent update.
  for (Index j = 0; j < m; j++) {
    Index k = (Index) I.elems[j];

    // Copy the kth column of matrix X.
    copyColumn(X.elems,x,k,n);

    // Perform the update.
    varbvsnormupdate(x,xy.elems[k],d.elems[k],sigma,sa,logodds.elems[k],
		     alpha.elems+k,mu.elems+k,Xr.elems,n);
  }

  // Free the dynamically allocated memory.
  free(x);
}
