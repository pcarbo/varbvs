varbvsmixupdate <- function (X, sigma, sa1, sa2, logodds, xy, d,
                             alpha0, mu10, mu20, Xr0, S) {
  # Runs a single iteration of the coordinate ascent updates to
  # maximize the variational lower bound for Bayesian variable
  # selection mixture model with linear regression. It adjusts the
  # fully-factorized variational approximation to the posterior
  # distribution of the coefficients in a linear regression model of a
  # continuous outcome (quantitiative trait), with normal mixture
  # priors on the coefficients.
  
  # CHECK THE INPUTS.
  # Check input X.
  if (!is.double(X) || !is.matrix(X))
    stop("Input argument 'X' must be a double-precision matrix")
  
  # Get the number of samples (n), the number of variables (p), and
  # the number of updates to execute (m).
  n <- nrow(X)
  p <- ncol(X)
  m <- length(S)

  # Check inputs sigma and sa.
  if (!is.scalar(sigma) || !is.scalar(sa1) || !is.scalar(sa2))
    stop("Input arguments 'sigma', 'sa1' and 'sa2' must be scalars")

  # Check input logodds.
  if (length(logodds) == 1)
    logodds <- rep(logodds,p)
  if (length(logodds) != p)
    stop("Input 'logodds' must be a scalar or a vector of length 'p'")

  # Check inputs xy and d.
  if (length(xy) != p || length(d) != p)
    stop("Inputs 'xy' and 'd' must be vectors of length 'p'")

  # Check inputs alpha0, mu10 and mu20.
  if (length(alpha0) != p || length(mu10) != p || length(mu20) != p)
    stop("Inputs 'alpha0', 'mu10' and 'mu20' must be vectors of length 'p'")

  # Check input Xr0.
  if (length(Xr0) != n)
    stop("Input 'Xr0' must be a vector of length 'n'")

  # Check input S.
  if (sum(S < 1 | S > p) > 0)
      stop("Input 'S' contains invalid variable indices")

  # Execute the C routine, and return the results in a list object.
  # The only components of the list that change are alpha, mu1, mu2
  # and Xr.  We need to subtract 1 from the indices because R vector
  # start at 1, but C arrays start at 0. Note that I do not attempt to
  # coerce X here; if X is large, it could use a lot of memory to
  # duplicate this matrix. For this same reason, I set DUP = FALSE so
  # that the input arguments are not duplicated.
  out <- .C("varbvsmixupdateR",
            n       = as.integer(n),      # Number of samples.
            m       = as.integer(m),      # Number of updates.
            X       = X,                  # Matrix of samples.
            sigma   = as.double(sigma),   # Residual variance.
            sa1     = as.double(sa1),     # Prior variance of coefficients.
            sa2     = as.double(sa2),     # Prior variance of coefficients.
            logodds = as.double(logodds), # Prior log-odds.
            xy      = as.double(xy),      # xy = X'*y.
            d       = as.double(d),       # d = diag(X'*X).
            alpha   = as.double(alpha0),  # Posterior inclusion prob's.
            mu1     = as.double(mu10),    # Posterior mean coefficients.
            mu2     = as.double(mu20),    # Posterior mean coefficients.
            Xr      = as.double(Xr0),     # Xr = X*(alpha*mu).
            S       = as.integer(S-1),    # Updates to perform.
            DUP     = FALSE)
  return(list(alpha = out$alpha,
              mu1   = out$mu1,
              mu2   = out$mu2,
              Xr    = out$Xr))
}
