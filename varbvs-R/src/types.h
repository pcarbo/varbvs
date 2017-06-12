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
#ifndef INCLUDE_TYPES
#define INCLUDE_TYPES

#ifdef MATLAB_MEX_FILE

// These type definitions are used to build the C routines for MATLAB.
#include "matrix.h"

typedef mwSize  Size;
typedef mwIndex Index;
typedef float   MatrixElem;

#else
#include <R.h>
#include <Rinternals.h>

// These type definitions are used to build the C routines for R. 
typedef R_xlen_t Size;
typedef R_xlen_t Index;
typedef double   MatrixElem;

#endif

#endif
