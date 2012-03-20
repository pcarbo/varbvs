#ifndef INCLUDE_VECTOROPS
#define INCLUDE_VECTOROPS

#include "types.h"

// Function declarations.
// -----------------------------------------------------------------
// Copy entries of one vector to another vector.
void copy (const double* source, double* dest, Size n);

// Get a pointer to column j of matrix X. The entries in a single
// column of the matrix are assumed to be stored consecutively in
// memory. Input n is the number of rows in the matrix.
const MatrixElem* getColumn (const MatrixElem* X, Index j, Size n);

// Copy column j of matrix X. The entries in a single column of the
// matrix are assumed to be stored consecutively in memory. Input n
// is the number of rows in the matrix.
void copyColumn (const MatrixElem* X, double* y, Index j, Size n);

// Set the entries of the vector to a.
void setVector (double* x, Size n, double a);

// Return the sum of the entries in the vector.
double sum (const double* x, Size n);

// Return the largest entry in the vector.
double max (const double* x, Size n);

// Add a to all the entries in vector x, and store the result in vector y.
void add (double* y, double a, const double* x, Size n);

// Return the dot product of two vectors.
double dot (const double* x, const double* y, Size n);

// Compute x'*D*y, the dot product of x and y scaled by D = diag(d).
double dotscaled (const double* x, const double* y, const double* d, Size n);

#endif
