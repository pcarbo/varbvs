# SUMMARY
# -------
# This file contains various functions to implement the variational
# methods for large-scale Bayesian variational selection. Here is an
# overview of the functions defined in this file:
#
#   tf2yn(x)
#   caterase(s)
#   dot(x,y)
#   norm2(x)
#   qnorm(x,a)
#   rep.col(x,n)
#   rep.row(x,n)
#   rand(m,n)
#   randn(m,n)
#   diagsq(X,a)
#   diagsqt(X,a)
#   diagsq2(X,A)
#   sigmoid(x)
#   logit(x)
#   logpexp(x)
#   logsigmoid(x)
#   slope(x)
#   int.gamma(logodds,alpha)
#   int.klbeta(alpha,mu,s,sa)
#   betavar(p,mu,s)
#   normalizelogweights(logw)
#
# Shorthand for machine precision.
eps <- .Machine$double.eps

# FUNCTION DEFINITIONS
# ----------------------------------------------------------------------
tf2yn <- function (x) {
  if (x)
    return("yes")
  else
    return("no")
}

# ----------------------------------------------------------------------
# Output the string using 'cat', then move the cursor back to the
# beginning of the string so that subsequent output will overwrite
# this string.
caterase <- function (s)
  cat(s,rep("\b",nchar(s)),sep = "")

# ----------------------------------------------------------------------
# Return the dot product of vectors x and y.
dot <- function (x,y)
  sum(x*y)

# ----------------------------------------------------------------------
# Return the quadratic norm (2-norm) of vector x.
norm2 <- function (x)
  sqrt(dot(x,x))

# ----------------------------------------------------------------------
# When input a is a matrix, this function returns the quadratic norm
# of vector x with respect to matrix a. When a is not a matrix, this
# function returns the norm of x with respect to A = diag(a). For a
# definition of the quadratic norm, see p. 635 of Convex Optimization
# (2004) by Boyd & Vandenberghe.
qnorm <- function (x, a) {
  if (is.matrix(a))
    y <- sqrt(c(x %*% a %*% x))
  else
    y <- sqrt(dot(x*a,x))
  return(y)
}

# ----------------------------------------------------------------------
# Replicate vector x to create an m x n matrix, where m = length(x).
rep.col <- function (x, n)
  matrix(x,length(x),n,byrow = FALSE)

# ----------------------------------------------------------------------
# Replicate vector x to create an n x m matrix, where m = length(x).
rep.row <- function (x, n)
  matrix(x,n,length(x),byrow = TRUE)

# ----------------------------------------------------------------------
# Return matrix containing pseudorandom values drawn from the standard
# uniform distribution.
rand <- function (m, n)
  matrix(runif(m*n),m,n)

# ----------------------------------------------------------------------
# Return matrix containing pseudorandom values drawn from the standard
# normal distribution.
randn <- function (m, n)
  matrix(rnorm(m*n),m,n)

# ----------------------------------------------------------------------
# diagsq(X) is the same as diag(X'*X), but the computation is done more
# efficiently, and without having to store an intermediate matrix of the
# same size as X. diag(X,a) efficiently computes diag(X'*diag(a)*X).
#
# This function calls "diag_Call", a function compiled from C code,
# using the .Call interface. To load the C function into R, first
# build the "shared object" (.so) file using the following command in
# the "src" directory: R CMD SHLIB diagsqr.c diagsq.c misc.c. Next,
# load the shared objects into R using the R function dyn.load:
# dyn.load("../src/diagsqr.so").
diagsq <- function (X, a = NULL) {

  # If input a is not provided, set it to a vector of ones.
  if (is.null(a))
    a <- rep(1,nrow(X))
  else
    a <- as.double(a)

  # Initialize the result.
  y <- rep(0,ncol(X))
  
  # Execute the C routine using the .Call interface. The main reason
  # for using .Call interface is that there is less of a constraint on
  # the size of the input matrices.
  out <- .Call("diagsq_Call",X = X,a = a,y = y)
  return(y)
}

# ----------------------------------------------------------------------
# diagsqt(X) returns the same result as diag(X*X'), but the
# computation is done more efficiently, and without having to store an
# intermediate result of the same size as X. diagsqt(X,a) efficiently
# computes diag(X*diag(a)*X').
#
# This function calls "diagt_Call", a function compiled from C code,
# using the .Call interface. See function "diagsq" for instructions on
# compiling and loading this C function into R.
diagsqt <- function (X, a = NULL) {

  # If input a is not provided, set it to a vector of ones.
  if (is.null(a))
    a <- rep(1,ncol(X))
  else
    a <- as.double(a)

  # Initialize the result.
  y <- rep(0,nrow(X))
  
  # Execute the C routine using the .Call interface. The main reason
  # for using .Call interface is that there is less of a constraint on
  # the size of the input matrices.
  out <- .Call("diagsqt_Call",X = X,a = a,y = y)
  return(y)
}

# ----------------------------------------------------------------------
# diagsq2(X) is the same as diag(X'*A*X), but the computation is done
# more efficiently.
diagsq2 <- function (X, A)
  rowSums((X %*% A) * X)

# ----------------------------------------------------------------------
# sigmoid(x) returns the sigmoid of the elements of x. The sigmoid
# function is also known as the logistic link function. It is the
# inverse of logit(x).
sigmoid <- function (x)
  1/(1 + exp(-x))

# ----------------------------------------------------------------------
# logit(x) returns the logit of the elements of X. It is the inverse of
# sigmoid(x).
logit <- function (x)
  log((x + eps)/((1 - x) + eps))
  
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
# slope(x) returns (sigmoid(x) - 1/2)/x, the slope of the conjugate to the
# log-sigmoid function at x, times 2. For details, see Bishop (2006), or the
# Bayesian Analysis paper. This is useful for working with the variational
# approximation for logistic regression.
slope <- function (x)
  (sigmoid(x) - 0.5)/(x + eps)

# ----------------------------------------------------------------------
# Computes an integral that appears in the variational lower bound of
# the marginal log-likelihood. This integral is the expectation on the
# prior inclusion probabilities taken with respect to the variational
# approximation. This returns the same result as sum(alpha*log(q) +
# (1-alpha)*log(1-q)).
int.gamma <- function (logodds, alpha)
  sum((alpha-1)*logodds + logsigmoid(logodds))

# ----------------------------------------------------------------------
# Computes an integral that appears in the variational lower bound of
# the marginal log-likelihood. This integral is the negative K-L
# divergence between the approximating distribution and the prior of
# the coefficients. Note that this sa is not the same as the sa used
# as an input to varbvsnorm.
int.klbeta <- function (alpha, mu, s, sa)
  (sum(alpha) + dot(alpha,log(s/sa)) - dot(alpha,s + mu^2)/sa)/2 -
    dot(alpha,log(alpha + eps)) - dot(1 - alpha,log(1 - alpha + eps))

# ----------------------------------------------------------------------
# Compute the variance of X, in which X is drawn from the normal
# distribution with probability p, and X is 0 with probability
# 1-p. Inputs mu and s specify the mean and variance of the normal
# density. Inputs p, mu and s must be arrays of the same
# dimension. This function is useful for calculating the variance of
# the coefficients under the fully-factorized variational
# approximation.
#
# Note that this is the same as 
# 
#    v = p*(s + mu^2) - (p*mu)^2.
#
betavar <- function (p, mu, s)
  p*(s + (1 - p)*mu^2)

# ----------------------------------------------------------------------
# normalizelogweights(logw) takes as input an array of unnormalized
# log-importance weights logw and returns normalized importance weights such
# that the sum of the normalized importance weights is equal to one.
normalizelogweights <- function (logw) {

  # Guard against underflow or overflow by adjusting the log-importance
  # weights so that the largest importance weight is one.
  c <- max(logw)
  w <- exp(logw - c)

  # Normalize the importance weights.
  return(w/sum(w))
}
