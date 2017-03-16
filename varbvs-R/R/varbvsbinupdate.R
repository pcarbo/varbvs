# Execute a single iteration of the coordinate ascent updates to
# maximize the variational lower bound for Bayesian variable selection
# in logistic regression.
#
# Input X is an n x p matrix of observations of the variables (or
# features), where n is the number of samples, and p is the number of
# variables. Input y contains samples of the binary outcome; it is a
# vector of length n.
#
# Input sa specifies the prior variance of the coefficients. Input
# logodds is the prior log-odds of inclusion for each variable. It
# must be a vector of length p. Note that a residual variance
# parameter (sigma) is not needed to model a binary outcome. See
# function 'updatestats_varbvsbin' for more information about input
# 'stats'.
#
# Inputs alpha0, mu0 are the current parameters of the variational
# approximation; under the variational approximation, the ith
# regression coefficient is normal with probability alpha0[i], and
# mu0[i] is the mean of the coefficient given that it is included in
# the model. Input Xr0 must be Xr0 = X*(alpha0*mu0).
#
# Input i specifies the order in which the coordinates are updated. It
# may be a vector of any length. Each entry of i must be an integer
# between 1 and p.
#
# There are three outputs. Output vectors alpha and mu are the updated
# variational parameters, and Xr = X*(alpha*mu). The computational
# complexity is O(n*length(i)).
#
# This function calls "varbvsbinupdate_Call", a function compiled from
# C code, using the .Call interface. See the comments accompanying
# function 'varbvsnormupdate' for instructions on building and loading
# the shared objects (.so) file into R.
varbvsbinupdate <- function (X, sa, logodds, stats, alpha0, mu0, Xr0, i) {

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
  out <- .Call(C_varbvsbinupdate_Call,X = X,sa = as.double(sa),
               logodds = as.double(logodds),d = as.double(stats$d),
               xdx = as.double(stats$xdx),xy = as.double(stats$xy),
               xd = as.double(stats$xd),alpha = alpha,mu = mu,Xr = Xr,
               i = as.integer(i-1))
  return(list(alpha = alpha,mu = mu,Xr = Xr))
}
