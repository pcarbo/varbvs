# Shorthand constant for machine precision.
eps <- .Machine$double.eps

is.scalar <- function (x) {
  # Returns TRUE if and only if x is a (numeric) scalar.
  return(length(x) == 1 && is.numeric(x))
}

is.odd <- function (x) {
  # Returns TRUE if and only if x is odd.
  return(x %% 2)
}

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

diagsq <- function (X, a = NULL) {
  # diagsq(X) returns diag(X'*X).
  # diagsq(X,a) returns diag(X'*diag(a)*X).
  if (is.null(a)) {
    n <- nrow(X)
    a <- rep(1,n)
  }

  # Compute y = (X.^2)'*a.
  a <- c(a)
  y <- c(a %*% X^2)  
  return(y)
}

diagsqt <- function (X, a = NULL) {
  # diagsqt(x) returns diag(X*X').  
  # diagsqt(X,a) returns diag(X*diag(a)*X').
  if (is.null(a)) {
    n <- ncol(X)
    a <- rep(1,n)
  }

  #  Compute y = X^2*a.
  a <- c(a)
  y <- c(X^2 %*% a)
  return(y)
}

repmat <- function (A,m,n) {
  # Does the same thing as repmat(A,m,n) in MATLAB.
  if (!is.matrix(A))
    stop("Invalid 'A' argument")
  return(kronecker(matrix(1,m,n),A))
}

center.columns <- function (X) {
  # Centers the columns of matrix X so that the entries in each column
  # of X add up to zero.
  if (!is.matrix(X))
    stop("Invalid 'X' argument")
  mu <- matrix(colMeans(X),1,ncol(X))
  X  <- X - repmat(mu,nrow(X),1)
  return(X)
}


relerr <- function (x1, x2) {
  # Returns the absolute relative error.
  if (length(x1) == 0 || length(x2) == 0)
    y <- 0
  else
    y <- abs(x1 - x2) / (abs(x1) + abs(x2) + eps)
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
  
logpexp <- function (x) {
  # logpexp(x) returns log(1 + exp(x)). For large x, logpexp(x) should
  # be approximately x. The computation is performed in a numerically
  # stable manner.

  # For large entries, log(1 + exp(x)) is effectively the same as x.
  y <- x

  # Find entries of x that are not large. For these entries, compute
  # log(1 + exp(x)).
  S    <- which(x < 16)
  y[S] <- log(1 + exp(x[S]))
  return(y)
}

logsigmoid <- function (x) {
  # Returns the logarithm of the sigmoid. Use this instead of
  # log(sigmoid(x)) to avoid loss of numerical precision.
  return(-logpexp(-x))
}

