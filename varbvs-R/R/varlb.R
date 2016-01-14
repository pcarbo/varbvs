intlinear <- function (Xr, d, y, sigma, alpha, mu, s) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the expectation
  # of the linear regression log-likelihood taken with respect to the
  # variational approximation. Inputs Xr and d must be equal to Xr =
  # X*(alpha*mu) and d = diag(X'*X). For a description of the
  # remaining inputs, see ‘varbvsoptimize’.
  n <- length(y)
  f <- (-n/2*log(2*pi*sigma) - norm2(y - Xr)^2/(2*sigma)
        - dot(d,betavar(alpha,mu,s))/(2*sigma))
  return(f)
}
