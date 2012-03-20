#ifndef INCLUDE_DIAGSQT
#define INCLUDE_DIAGSQT

#include "types.h"

// Function declarations.
// -----------------------------------------------------------------
// Compute X.^2*a and store the result in vector y. X is an m x n
// matrix in which the entries in each column are stored consecutively
// in memory.
void diagsqt (const MatrixElem* X, const double* a, double* y, 
	      Size m, Size n);

#endif
