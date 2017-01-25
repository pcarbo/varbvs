#include <Rcpp.h>
using namespace Rcpp;

// Return the sigmoid function at x.
// [[Rcpp::export]]
double sigmoidcpp (double x) {
  return(1/(1 + exp(-x)));
}

//' Variational update of normal varbvs.
//'
//' @param X An n by p matrix of predictors. The rows index the observations and the columns index the variables.
//' @param sigma A positive numeric. The residual variance.
//' @param sa A positive numeric. The prior variance.
//' @param logodds A p-vector that contains the prior log-odds of inclusion for each variable.
//' @param xy A p-vector containing \code{t(X) \%*\% y}, where \code{y} is a response vector.
//' @param d The diagonal of the Gram matrix, \code{d = diag(t(X) \%*\% X)}
//' @param alpha A p-vector of inclusion probabilities for each variable.
//' @param mu A p-vector of means for each coefficient.
//' @param Xr An n-vector with elements \code{X \%*\% (alpha * mu)}
//' @param i A p-vector of integers that give the order of updating the coefficients.
//'
//' @export
//'
//' @author David Gerard
//'
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
    alpha[j] = sigmoidcpp(logodds[j] + (log(s/(sa*sigma)) + pow(mu[j],2)/s)/2);

    // Update Xr = X*r.
    Xr = Xr + (alpha[j] * mu[j] - r) * X.column(j);
  }
}
