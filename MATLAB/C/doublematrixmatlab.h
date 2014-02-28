#ifndef INCLUDE_DOUBLEMATRIXMATLAB
#define INCLUDE_DOUBLEMATRIXMATLAB

#include "types.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Type definitions.
// -----------------------------------------------------------------
// A dense matrix with single-precision floating-point entries.
typedef struct {
  Size    nr;     // Number of rows.
  Size    nc;     // Number of columns.
  double* elems;  // Entries of matrix.
} DoubleMatrix;

// Function declarations.
// -----------------------------------------------------------------
// Get a double-precision floating point matrix from a MATLAB array. 
DoubleMatrix getDoubleMatrix (const mxArray* ptr);

#endif
