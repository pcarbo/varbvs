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
