#ifndef INCLUDE_DOUBLEMATRIXMEX
#define INCLUDE_DOUBLEMATRIXMEX

#include "types.h"
#include "mex.h"
#include "matrix.h"

// TYPE DEFINITIONS
// -----------------------------------------------------------------
// A dense matrix with single precision entries.
typedef struct {
  Size    nr;     // Number of rows.
  Size    nc;     // Number of columns.
  double* elems;  // Entries of matrix.
} DoubleMatrix;

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Get a double precision floating point matrix from a MATLAB array. 
DoubleMatrix getDoubleMatrix (const mxArray* ptr);

// Get a pointer to column j of matrix X. The entries in a single
// column of the matrix are assumed to be stored consecutively in
// memory. Input n is the number of rows in the matrix.
const double* getDoubleColumn (const double* X, Index j, Size n);

// Copy column j of matrix X. The entries in a single column of the
// matrix are assumed to be stored consecutively in memory. Input n
// is the number of rows in the matrix.
void copyDoubleColumn (const double* X, double* y, Index j, Size n);

#endif
