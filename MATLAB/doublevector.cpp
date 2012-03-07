#include "doublevector.h"
#include "common.h"
#include <string.h>

// Function definitions.
// -----------------------------------------------------------------
// Allocate memory for a double-precision vector of length n. 
doublevector newdoublevector (mwSize n) {
  doublevector result;  // The return value.
  result.n     = n;
  result.elems = new double[n];
  return result;
}

// Free the dynamically allocated memory associated with a
// doublevector structure.
void free (doublevector& x) {
  delete[] x.elems;
}

// Copy the vector elements from the source to the destination. It is
// assumed that sufficient memory has been allocated for the
// destination vector.
void copy (const doublevector& source, doublevector& dest) {
  dest.n = source.n;
  memcpy(dest.elems,source.elems,sizeof(double)*source.n);
}

// Get a double-precision vector from a MATLAB array.
doublevector getdoublevector (const mxArray* ptr) {
  doublevector result;  // The return value.
  result.n     = mxGetNumberOfElements(ptr);
  result.elems = mxGetPr(ptr);
  return result;
}

// Create a MATLAB array and return the corresponding column vector.
doublevector creatematlabvector (mwSize n, mxArray*& ptr) {
  ptr = mxCreateDoubleMatrix(n,1,mxREAL);
  return getdoublevector(ptr);
}

// Set the entries of vector x to a.
void setvector (doublevector& x, double a) {
  for (mwIndex i = 0; i < x.n; i++)
    x.elems[i] = a;
}

// Return the sum of the entries of vector x.
double sum (const doublevector& x) {
  return sum(x.elems,x.n);
}

// Compute y = y + a*x.
void add (doublevector& y, double a, const doublevector& x) {
  for (mwIndex i = 0; i < x.n; i++)
    y.elems[i] += a * x.elems[i];
}

// Returns dot product of vectors x and y.
double dot (const doublevector& x, const doublevector& y) {
  double z = 0;  // The return value.
  for (mwIndex i = 0; i < x.n; i++)
    z += x.elems[i] * y.elems[i];
  return z;
}

// Returns x'*D*y, the dot product of x and y scaled by matrix D, where 
// D = diag(d).
double dot (const doublevector& x, const doublevector& y, 
	    const doublevector& d) {
  double z = 0;  // The return value.
  for (mwIndex i = 0; i < x.n; i++)
    z += x.elems[i] * y.elems[i] * d.elems[i];
  return z;
}
