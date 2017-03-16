# Executes a single iteration of the coordinate ascent updates to
# maximize the variational lower bound for Bayesian variable selection
# in logistic regression, allowing for additional covariates. See
# varbvsbinupdate for more details.
varbvsbinzupdate <- function (X, sa, logodds, stats, alpha0, mu0, Xr0, i) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  # Check input X.
  if (!is.double(X) || !is.matrix(X))
    stop("Input argument 'X' must be a double-precision matrix")

  # Check input sa.
  if (length(sa) != 1)
    stop("Input sa must be a scalar")

  # Check input logodds, alpha0 and mu0.
  if (!(length(logodds) == p & length(alpha0) == p & length(mu0) == p))
    stop("logodds, alpha0 and mu0 must have length = ncol(X)")

  # Check input Xr0.
  if (length(Xr0) != n)
    stop("length(Xr0) must be equal to nrow(X)")

  # Check input i.
  if (sum(i < 1 | i > p) > 0)
    stop("Input i contains invalid variable indices")

  # Initialize storage for the results.
  alpha <- c(alpha0)
  mu    <- c(mu0)
  Xr    <- c(Xr0)

  # Execute the C routine using the .Call interface, and return the
  # updated variational parameters statistics in a list object. The
  # main reason for using the .Call interface is that there is less of
  # a constraint on the size of the input matrices. The only
  # components that change are alpha, mu and Xr. Note that I need to
  # subtract 1 from the indices because R vectors start at 1, and C
  # arrays start at 0.
  out <- .Call(C_varbvsbinzupdate_Call,X = X,sa = as.double(sa),
               logodds = as.double(logodds),d = as.double(stats$d),
               xdx = as.double(stats$xdx),xy = as.double(stats$xy),
               dzr = stats$dzr,alpha = alpha,mu = mu,Xr = Xr,
               i = as.integer(i-1))
  return(list(alpha = alpha,mu = mu,Xr = Xr))
}
