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
