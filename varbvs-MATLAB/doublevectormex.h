#ifndef INCLUDE_DOUBLEVECTORMEX
#define INCLUDE_DOUBLEVECTORMEX

#include "types.h"

// These include files have a bunch of definitions to interface C
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// TYPE DEFINITIONS
// -----------------------------------------------------------------
// A vector with double precision entries.
typedef struct {
  Size    n;      // Size of vector.
  double* elems;  // Vector entries.
} DoubleVector;

// FUNCTION DECLARATIONS
// -----------------------------------------------------------------
// Get a double precision vector from a MATLAB array.
DoubleVector getDoubleVector (const mxArray* ptr);

// Create a column vector in MATLAB array.
DoubleVector createMatlabVector (Size n, mxArray** ptr);

// Copy entries of one vector to another vector.
void copyDoubleVector (const DoubleVector source, DoubleVector dest);

#endif
