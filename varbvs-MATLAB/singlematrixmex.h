#ifndef INCLUDE_SINGLEMATRIXMEX
#define INCLUDE_SINGLEMATRIXMEX

#include "types.h"
#include "mex.h"
#include "matrix.h"

// TYPE DEFINIITIONS
// -----------------------------------------------------------------
// A dense matrix with single precision floating point entries.
typedef struct {
  Size        nr;     // Number of rows.
  Size        nc;     // Number of columns.
  MatrixElem* elems;  // Entries of matrix.
} SingleMatrix;

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Get a single precision matrix from a MATLAB array. 
SingleMatrix getSingleMatrix (const mxArray* ptr);

#endif
