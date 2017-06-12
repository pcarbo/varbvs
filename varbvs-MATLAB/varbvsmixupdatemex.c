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
// For a description of this C code, see varbvsmixupdate.m.
#include "types.h"
#include "misc.h"
#include "doublevectormex.h"
#include "singlematrixmex.h"
#include "doublematrixmex.h"
#include "varbvs.h"
#include "mex.h"
#include "matrix.h"

void mexFunction (int nlhs, mxArray* plhs[], 
		  int nrhs, const mxArray* prhs[]) {

  // (1) GET INPUTS
  // --------------
  const SingleMatrix X      = getSingleMatrix(prhs[0]);
  const double       sigma  = *mxGetPr(prhs[1]);
  const DoubleVector sa     = getDoubleVector(prhs[2]);
  const DoubleVector q      = getDoubleVector(prhs[3]);
  const DoubleVector xy     = getDoubleVector(prhs[4]);
  const DoubleVector d      = getDoubleVector(prhs[5]);
  const DoubleMatrix alpha0 = getDoubleMatrix(prhs[6]);
  const DoubleMatrix mu0    = getDoubleMatrix(prhs[7]);
  const DoubleVector Xr0    = getDoubleVector(prhs[8]);
  const DoubleVector I      = getDoubleVector(prhs[9]);
  const double       eps    = *mxGetPr(prhs[10]);
  
  // Get the number of samples (n), the number of variables (p), the
  // number of mixture components (K), and the number of coordinate
  // ascent updates (m).
  const Size n = X.nr;
  const Size p = X.nc;
  const Size k = q.n;
  const Size m = I.n;

  // (2) INITIALIZE OUTPUTS
  // ----------------------
  DoubleMatrix alpha = createMatlabDoubleMatrix(k,p,&plhs[0]);
  DoubleMatrix mu    = createMatlabDoubleMatrix(k,p,&plhs[1]);
  DoubleVector Xr    = createMatlabVector(n,&plhs[2]);

  copyDoubleMatrix(alpha0,alpha);
  copyDoubleMatrix(mu0,mu);
  copyDoubleVector(Xr0,Xr);

  // These variables provide storage for a column of matrix X, and
  // other intermediate calculations for the co-ordinate ascent
  // updates below.
  double* x    = malloc(sizeof(double)*n);
  double* s    = malloc(sizeof(double)*k);
  double* logw = malloc(sizeof(double)*k);

  // (3) RUN COORDINATE ASCENT UPDATES
  // ---------------------------------
  // Repeat for each coordinate ascent update.
  for (Index j = 0; j < m; j++) {
    Index i = (Index) I.elems[j];

    // Copy the ith column of matrix X.
    copyColumn(X.elems,x,i,n);

    // Perform the update.
    varbvsmixupdate(x,xy.elems[i],d.elems[i],sigma,sa.elems,q.elems,
		    getDoubleColumn(alpha.elems,i,k),
		    getDoubleColumn(mu.elems,i,k),
		    Xr.elems,s,logw,n,k,eps);
  }

  // Free the dynamically allocated memory.
  free(x);
  free(s);
  free(logw);
}
