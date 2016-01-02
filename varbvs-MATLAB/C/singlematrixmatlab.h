#ifndef INCLUDE_SINGLEMATRIXMATLAB
#define INCLUDE_SINGLEMATRIXMATLAB

#include "types.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Type definitions.
// -----------------------------------------------------------------
// A dense matrix with single-precision floating-point entries.
typedef struct {
  Size        nr;     // Number of rows.
  Size        nc;     // Number of columns.
  MatrixElem* elems;  // Entries of matrix.
} SingleMatrix;

// Function declarations.
// -----------------------------------------------------------------
// Get a single precision floating point matrix from a MATLAB array. 
SingleMatrix getSingleMatrix (const mxArray* ptr);

#endif
