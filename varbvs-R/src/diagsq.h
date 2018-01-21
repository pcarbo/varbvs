// Part of the varbvs package, https://github.com/pcarbo/varbvs
//
// Copyright (C) 2012-2018, Peter Carbonetto
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
#ifndef INCLUDE_DIAGSQ
#define INCLUDE_DIAGSQ

#include "types.h"

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Compute (X.^2)'*a and store the result in vector y. X is an m x n
// matrix in which the entries in each column are stored consecutively
// in memory.
void diagsq (const MatrixElem* X, const double* a, double* y, 
	     Size m, Size n);

// Compute X.^2*a and store the result in vector y. X is an m x n
// matrix in which the entries in each column are stored consecutively
// in memory.
void diagsqt (const MatrixElem* X, const double* a, double* y, Size m, Size n);

#endif
