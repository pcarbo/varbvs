#ifndef INCLUDE_TYPES
#define INCLUDE_TYPES

#ifdef MATLAB_MEX_FILE

// This include file has a bunch of definitions to interface C
// routines to MATLAB.
#include "matrix.h"

typedef mwSize  Size;
typedef mwIndex Index;
typedef float   MatrixElem;

#endif

#endif
