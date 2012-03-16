#ifndef INCLUDE_DOUBLEVECTOR
#define INCLUDE_DOUBLEVECTOR

// These include files have a bunch of definitions to interface C++
// routines to MATLAB.
#include "mex.h"
#include "matrix.h"

// Type definitions.
// -----------------------------------------------------------------
// A vector in double precision format.
typedef struct DoubleVector {
  mwSize  n;      // Length of vector.
  double* elems;  // Vector elements.
} doublevector;

// Function declarations.
// -----------------------------------------------------------------
// Allocate memory for a double-precision vector of length n. 
doublevector newdoublevector (mwSize n);

// Free the dynamically allocated memory associated with a
// doublevector structure.
void free (doublevector& x);

// Copy the vector elements from the source to the destination. It is
// assumed that sufficient memory has been allocated for the
// destination vector.
void copy (const doublevector& source, doublevector& dest);

// Get a double-precision vector from a MATLAB array.
doublevector getdoublevector (const mxArray* ptr);

// Create a MATLAB array and return the corresponding column vector.
doublevector creatematlabvector (mwSize n, mxArray*& ptr);

// Set the entries of vector x to a.
void setvector (doublevector& x, double a);

// Return the sum of the entries of vector x.
double sum (const doublevector& x);

// Compute y = y + a*x.
void add (doublevector& y, double a, const doublevector& x);

// Returns dot product of vectors x and y.
double dot (const doublevector& x, const doublevector& y);

// Returns x'*D*y, the dot product of x and y scaled by matrix D, where 
// D = diag(d).
double dot (const doublevector& x, const doublevector& y, 
	    const doublevector& d);

#endif
