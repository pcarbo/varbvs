# Execute a single iteration of the coordinate ascent updates to
# maximize the variational lower bound for Bayesian variable selection
# in linear regression.
#
# Input X is an n x p matrix of observations of the variables (or
# features), where n is the number of samples, and p is the number of
# variables. Input xy = X'*y, where y is the vector of samples of the
# continuous outcome.
#
# Inputs sigma, sa and logodds specify other model parameters. sigma
# and sa are scalars. sigma specifies the variance of the residual,
# and sa is the prior variance of the regression coefficients (scaled
# by sigma). Input logodds is the prior log-odds of inclusion for each
# variable. It must be a vector of length p.
#
# Inputs alpha0, mu0 are the current parameters of the variational
# approximation. Under the variational approximation, the ith
# regression coefficient is included in the model with probability
# alpha0[i], and mu0(i) is the mean of the coefficient given that it
# is included in the model. Inputs Xr0 and d must be Xr0 =
# X*(alpha0*mu0) and d = diag(X'*X).
#
# Input i specifies the order in which the coordinates are updated. It
# may be a vector of any length. Each entry of i must be an integer
# between 1 and p.
#
# There are 3 outputs. Output vectors alpha and mu are the updated
# variational parameters, and Xr = X*(alpha*mu). The computational
# complexity is O(n*length(i)).
#

# When algorithm.version = ".Call", this function calls "varbvsnormupdate_Call",
# a function compiled from C code, using the .Call interface. To load
# the C function into R, first build the "shared object" (.so) file
# using the following command in the "src" directory: R CMD SHLIB
# varbvsr.c varbvs.c misc.c. Next, load the shared objects into R
# using the R function dyn.load: dyn.load("../src/varbvsr.so").

#
# Add instructions for updating the 
    #
    # TO DO: Add instructions for updating the Rcpp code.
    #

varbvsnormupdate <-
    function (X, sigma, sa, logodds, xy, d, alpha0, mu0, Xr0, i,
              algorithm.version = c(".Call","Rcpp","R")) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  # Specify the algorithm implementation.
  algorithm.version <- match.arg(algorithm.version)

  # Check input X.
  if (!is.double(X) || !is.matrix(X))
    stop("Input X must be a double-precision matrix")

  # Check inputs sigma and sa.
  if (length(sigma) != 1 | length(sa) != 1)
    stop("Inputs sigma and sa must be scalars")

  # Check input logodds, xy, d, alpha0 and mu0.
  if (!(length(logodds) == p & length(xy) == p & length(d) == p &
       length(alpha0) == p & length(mu0) == p))
    stop("logodds, xy, d, alpha0 and mu0 must have length = ncol(X)")

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

  if (algorithm.version == ".Call") {

    # Execute the C routine using the .Call interface, and return the
    # updated variational parameters statistics in a list object. The
    # main reason for using the .Call interface is that there is less of
    # a constraint on the size of the input matrices. The only
    # components that change are alpha, mu and Xr. Note that I need to
    # subtract 1 from the indices because R vectors start at 1, and C
    # arrays start at 0.
    out <- .Call("varbvsnormupdate_Call",X = X,sigma = as.double(sigma),
                 sa = as.double(sa),logodds = as.double(logodds),
                 xy = as.double(xy),d = as.double(d),alpha = alpha,mu = mu,
                 Xr = Xr,i = as.integer(i-1))
  } else if (algorithm.version == "Rcpp") {

    # Execute the C routine using the Rcpp interface.
    varbvsnormupdate_Rcpp(X = X,sigma = sigma,sa = sa,logodds = logodds,
                          xy = xy,d = d,alpha = alpha,mu = mu,Xr = Xr,
                          i = i-1)
  } else if (algorithm.version == "R") {

    # Repeat for each co-ordinate to update.
    for (j in i) {

      # Compute the variational estimate of the posterior variance.
      s <- sa*sigma/(sa*d[j] + 1)

      # Update the variational estimate of the posterior mean.
      r     <- alpha[j] * mu[j]
      mu[j] <- s/sigma * (xy[j] + d[j]*r - sum(X[,j]*Xr))

      # Update the variational estimate of the posterior inclusion
      # probability.
      alpha[j] <- sigmoid(logodds[j] + (log(s/(sa*sigma)) + mu[j]^2/s)/2)

      # Update Xr = X*r.
      Xr <- Xr + (alpha[j]*mu[j] - r) * X[,j]
    }
  } else
    stop("Invalid input argument 'algorithm.version' in varbvsnormupdate.")
  return(list(alpha = alpha,mu = mu,Xr = Xr))
}
