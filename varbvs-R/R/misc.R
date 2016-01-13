# SUMMARY
# -------

# This file contains various functions to implement the variational
# methods for large-scale Bayesian variational selection. Here is an
# overview of the functions defined in this file:
#
#   tf2yn(x)
#   diagsq(X,a)
#   logpexp(x)
#   logsigmoid(x)
#   int.gamma(logodds,alpha)
#
# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
tf2yn <- function (x) {
  if (x)
    return("yes")
  else
    return("no")
}

# ----------------------------------------------------------------------
# diagsq(X) is the same as diag(X'*X), but the computation is done more
# efficiently, and without having to store an intermediate matrix of the
# same size as X. diag(X,a) efficiently computes diag(X'*diag(a)*X).
diagsq <- function (X, a = NULL) {

  # If input a is not provided, set it to a vector of ones.
  if (is.null(a))
    a <- rep(1,nrow(X))
  else
    a <- c(a)
  
  # Compute the result using the efficient C routine.
  #
  # TO DO: Implement the efficient C routine.
  #
  return(a %*% X^2)
}

# ----------------------------------------------------------------------
# logpexp(x) returns log(1 + exp(x)). The computation is performed in a
# numerically stable manner. For large entries of x, log(1 + exp(x)) is
# effectively the same as x.
logpexp <- function (x) {
  y    <- x
  i    <- which(x < 16)
  y[i] <- log(1 + exp(x[i]))
  return(y)
}

# ----------------------------------------------------------------------
# Use this instead of log(sigmoid(x)) to avoid loss of numerical precision.
logsigmoid <- function (x)
  -logpexp(-x)

# ----------------------------------------------------------------------
# Computes an integral that appears in the variational lower bound of
# the marginal log-likelihood. This integral is the expectation on the
# prior inclusion probabilities taken with respect to the variational
# approximation. This returns the same result as sum(alpha*log(q) +
# (1-alpha)*log(1-q)).
int.gamma <- function (logodds, alpha)
  sum((alpha-1)*logodds + logsigmoid(logodds))

# ****** OLD STUFF ******

# Shorthand constant for machine precision.
eps <- .Machine$double.eps

dot <- function (x,y) {
  # Returns the dot product of vectors x and y.
  return(sum(x*y))
}

norm2 <- function (x) {
  # Returns the quadratic norm (2-norm) of vector x.
  return(sqrt(dot(x,x)))
}

qnorm <- function (x, a) {
  # Returns the quadratic norm of vector x with respect to positive
  # definite matrix with diagonal entries a. For a definition of the
  # quadratic norm, see p. 635 of Convex Optimization (2004) by Boyd &
  # Vandenberghe.
  return(sqrt(dot(x*a,x)))
}

diagsqt <- function (X, a = NULL) {
  # diagsqt(x) returns diag(X*X').  
  # diagsqt(X,a) returns diag(X*diag(a)*X').
  if (is.null(a)) {
    n <- ncol(X)
    a <- rep(1,n)
  }

  # Compute y = X^2*a.
  a <- c(a)
  y <- c(X^2 %*% a)
  return(y)
}

normalizelogweights <- function (logw) {
  # Takes as input an array of unnormalized log-importance weights and
  # returns normalized importance weights such that the sum of the
  # normalized importance weights is equal to one.

  # We guard against underflow or overflow by adjusting the log-importance
  # weights so that the largest importance weight is one.
  c <- max(logw)
  w <- exp(logw - c)

  # Normalize the importance weights.
  w <- w / sum(w)
}

betavar <- function (p, mu, s) {
  # betavar(p,mu,s) returns the variance of random vector X, in which
  # X[i] is drawn from the normal distribution with probability p[i],
  # and X[i] is zero with probability 1-p[i]. Inputs mu and s specify
  # the mean and variance of the normal density. Inputs p, mu and s
  # must be vectors of the same length. This function is useful for
  # calculating the variance of the coefficients under the
  # fully-factorized variational approximation.

  # Note that this is the same as 
  # 
  #    v = p*(s + mu^2) - (p*mu)^2.
  #
  return(p*(s + (1 - p)*mu^2))
}

sigmoid <- function (x) {
  # Returns the sigmoid of x. The sigmoid function is also known as
  # the logistic link function. It is the inverse of logit(x).
  return(1/(1 + exp(-x)))
}

logit <- function (x) {
  # The logit function, which is the reverse of the sigmoid function.
  return(log((x + eps)/((1 - x) + eps)))
}
  
