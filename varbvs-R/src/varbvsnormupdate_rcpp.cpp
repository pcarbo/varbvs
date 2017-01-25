#include <Rcpp.h>
using namespace Rcpp;

// Return the sigmoid function at x.
// [[Rcpp::export]]
double sigmoid_rcpp (double x) {
  return(1/(1 + exp(-x)));
}

// Execute a single coordinate ascent update to maximize the variational 
// lower bound for Bayesian variable selection in linear regression.
// [[Rcpp::export]]
void varbvsnormupdate_rcpp (const NumericMatrix& X, double sigma, double sa,
			    const NumericVector& logodds,
			    const NumericVector& xy, const NumericVector& d,
			    NumericVector& alpha, NumericVector& mu,
			    NumericVector& Xr, const IntegerVector& i) {

  // Cycle through the coordinate ascent updates.
  for (int iter = 0; iter < i.size(); iter++) {

    // j is the coordinate we update in this iteration.
    int j = i[iter];

    // Compute the variational estimate of the posterior variance.
    double s = sa*sigma/(sa*d[j] + 1);

    // Update the variational estimate of the posterior mean.
    double r = alpha[j] * mu[j];
    mu[j] = s/sigma*(xy[j] + d[j]*r - sum(X.column(j) * Xr));

    // Update the variational estimate of the posterior inclusion probability.
    alpha[j] = sigmoid_rcpp(logodds[j] + (log(s/(sa*sigma)) +
					  mu[j]*mu[j]/s)/2);

    // Update Xr = X*r.
    Xr = Xr + (alpha[j] * mu[j] - r) * X.column(j);
  }
}
