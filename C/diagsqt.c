#include "diagsqt.h"

// Function definitions.
// -----------------------------------------------------------------
// Compute X.^2*a and store the result in vector y.
void diagsqt (const MatrixElem* X, const double* a, double* y, 
	     Size m, Size n) {
  double t;  // An intermediate result.

  // Zero the entries of the result vector.
  setvector(y,0);

  // Repeat for each column of X.
  for (Index j = 0; j < n; j++, a++)

    // Repeat for each row of X.
    const double* yi = y;
    for (Index i = 0; i < m; i++, X++, yi++) {
      
      // Add X(i,j)^2 * a(j) to the ith entry of y.
      t    = (double) *X;
      *yi += t * t * (*a);
    }
}
