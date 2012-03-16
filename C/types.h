#ifndef INCLUDE_TYPES
#define INCLUDE_TYPES

#ifdef MATLAB_MEX_FILE

// This include file has a bunch of definitions to interface C
// routines to MATLAB.
#include "matrix.h"

// These definitions are used to build the C routines for MATLAB.
typedef mwSize  Size;
typedef mwIndex Index;
typedef float   MatrixElem;

#else

// These definitions are used to build the C routines for R.
typedef int    Size;
typedef int    Index;
typedef double MatrixElem;

#endif
