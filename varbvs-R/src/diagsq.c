#include "diagsq.h"
#include "misc.h"

// FUNCTION DEFINITIONS
// ---------------------------------------------------------------------
// Compute (X.^2)'*a and store the result in vector y.
void diagsq (const MatrixElem* X, const double* a, double* y, 
	     Size m, Size n) {
  double t;

  // Repeat for each column of X.
  for (Index j = 0; j < n; j++, y++) {
    
    // Initialize the jth entry of the result.
    *y = 0;

    // Repeat for each row of X.
    const double* ai = a;
    for (Index i = 0; i < m; i++, X++, ai++) {

      // Add X(i,j)^2 * a(i) to the jth entry of y.
      t   = (double) *X;
      *y += t * t * (*ai);
    }
  }
}

// ---------------------------------------------------------------------
// Compute X.^2*a and store the result in vector y.
void diagsqt (const MatrixElem* X, const double* a, double* y, Size m,
	      Size n) {
  double t;  // An intermediate result.

  // Zero the entries of the result vector.
  setVector(y,m,0);

  // Repeat for each column of X.
  for (Index j = 0; j < n; j++, a++) {

    // Repeat for each row of X.
    double* yi = y;
    for (Index i = 0; i < m; i++, X++, yi++) {
      
      // Add X(i,j)^2 * a(j) to the ith entry of y.
      t    = (double) *X;
      *yi += t * t * (*a);
    }
  }
}
