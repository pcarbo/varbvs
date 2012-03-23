# Shorthand constant for machine precision.
eps <- .Machine$double.eps

# NORMALIZELOGWEIGHTS(LOGW) takes as input an array of unnormalized
# log-importance weights LOGW and returns normalized importance weights such
# that the sum of the normalized importance weights is equal to one.
normalizelogweights <- function (logw) {

  # We guard against underflow or overflow by adjusting the log-importance
  # weights so that the largest importance weight is one.
  c <- max(c(logw))
  w <- exp(logw - c)

  # Normalize the importance weights.
  w <- w / sum(c(w))
}
  
# LOGINVGAMMA(X,A,B) returns the logarithm of the probability density
# function of the inverse gamma distribution at elements of X, with
# shape A and scale B.
loginvgamma <- function (x, a, b) {
  f <- a*log(b) - lgamma(a) - (a+1)*log(x+eps) - b/(x+eps)
  return(f)
}
  
# LOGPVE(A,X) returns the logarithm of the density function for X,
# given that R = A*X/(A*X + 1) is uniform on [0,1]. Input A must be a
# positive scalar. This is useful for calculating the prior
# probability of the variance parameter X, when we have a uniform
# prior on R, the proportion of variance explained.
logpve <- function (a, x) {
  r <- a*x/(a*x+1)          # Proportion of variance explained.
  f <- 2*log(r/x) - log(a)  # Log-density.
  return(f)
}
  
intlinear <- function (Xr, d, y, sigma, alpha, mu, s) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the expectation
  # of the linear regression log-likelihood taken with respect to the
  # variational approximation. Inputs Xr and d must be equal to Xr =
  # X*(alpha*mu) and d = diag(X'*X). For a description of the
  # remaining inputs, see function ‘varbvsoptimize’.
  n <- length(y)
  f <- -n/2*log(2*pi*sigma) - norm2(y - Xr)^2/(2*sigma) -
        dot(d,betavar(alpha,mu,s))/(2*sigma)
  return(f)
}
                      
intklbeta <- function (alpha, mu, s, sa) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the negative
  # Kullback-Leibler divergence between the approximating distribution
  # and the prior of the coefficients. Input ‘sa’ specifies the prior
  # variance of the coefficients. (It is not the same as the ‘sa’ used
  # as input to ‘varbvsoptimize’.) See function ‘varbvsoptimize’ for
  # details on the inputs to this function.
  f <- (sum(alpha) + dot(alpha,log(s/sa)) - dot(alpha,s + mu^2)/sa)/2 -
        dot(alpha,log(alpha + eps)) -
        dot(1 - alpha,log(1 - alpha + eps))
  return(f)
}

intgamma <- function (logodds, alpha) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the expectation
  # on the prior inclusion probabilities taken with respect to the
  # variational approximation.
  #
  # Args:
  #   logodds  Scalar, or a vector, specifying the prior log-odds.
  #   alpha    Mixture weights for the variational approximation.
  #
  # See function ‘varbvsoptimize’ for details on the input arguments.

  # This is the same as 
  #
  #    sum(alpha*log(q) + (1-alpha)*log(1-q)).
  #  
  return(sum((alpha-1) * logodds + logsigmoid(logodds)))
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

dot <- function (x,y) {
  # Return the dot product of vectors x and y.
  return(sum(x*y))
}

norm2 <- function (x) {
  # Return the quadratic norm (2-norm) of vector x.
  return(sqrt(dot(x,x)))
}

logsigmoid <- function (x) {
  # Returns the logarithm of the sigmoid. Use this instead of
  # log(sigmoid(x)) to avoid loss of numerical precision.
  return(-logpexp(-x))
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

logit <- function (x) {
  # The logit function, which is the reverse of the sigmoid function.
  return(log((x + eps)/((1 - x) + eps)))
}
  
diagsq <- function (X, a = NULL) {
  # diagsq(X) returns diag(X'*X).
  # diagsq(X,a) returns diag(X'*diag(a)*X).
  if (is.null(a)) {
    n <- nrow(X)
    a <- rep(1,n)
  }
  y <- c(a %*% X^2)  # y = X^2*a.
  return(y)
}

relerr <- function (x1, x2) {
  # Returns the absolute relative error.
  if (length(x1) == 0 || length(x2) == 0)
    y <- 0
  else
    y <- abs(x1 - x2) / (abs(x1) + abs(x2) + eps)
  return(y)
}

repmat <- function (A,m,n) {
  # Does the same thing as REPMAT(A,m,n) in MATLAB.
  return(kronecker(matrix(1,m,n),A))
}

center.columns <- function (X) {
  # Centers the columns of X so that the entries in each column of X
  # add up to zero.
  mu <- matrix(colMeans(X),1,ncol(X))
  mu <- repmat(mu,nrow(X),1)
  X  <- X - mu
  return(X)
}

grid3d <- function (x,y,z) {
  # Does the same thing as NDGRID(x,y,z) in MATLAB.

  # Get the number of entries in each of the inputs.
  nx <- length(x)
  ny <- length(y)
  nz <- length(z)

  # Initialize the 3-d outputs.
  X <- array(dim=c(nx,ny,nz))
  Y <- X
  Z <- X

  # Set the entries of the 3-d arrays.
  for (i in 1:nx)
    for (j in 1:ny)
      for (k in 1:nz) {
        X[i,j,k] <- x[i]
        Y[i,j,k] <- y[j]
        Z[i,j,k] <- z[k]
      }
  
  grid <- list(X = X,Y = Y,Z = Z)
  return(grid)
}

is.scalar <- function (x) {
  # Returns TRUE if and only if x is a scalar.
  return(length(x) == 1)
}

is.odd <- function (x) {
  # Returns TRUE if and only if x is odd.
  return(x %% 2)
}

