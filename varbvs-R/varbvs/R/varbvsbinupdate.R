varbvsbinupdate <- function (X, sa, logodds, stats, alpha0, mu0, Xr0, S) {
  # Runs a single iteration of the coordinate ascent updates to
  # maximize the variational lower bound for Bayesian variable
  # selection in logistic regression. It adjusts the fully-factorized
  # variational approximation to the posterior distribution of the
  # coefficients in a logistic regression model of a binary outcome or
  # trait, with spike and slab priors on the coefficients.

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
  if (!is.scalar(sa))
    stop("Input arguments 'sa' must be a scalar")

  # Check input logodds.
  if (length(logodds) == 1)
    logodds <- rep(logodds,p)
  if (length(logodds) != p)
    stop("Input 'logodds' must be a scalar or a vector of length 'p'")

  # Check input stats.
  if (!is.list(stats) || is.null(stats$u) || is.null(stats$xy) ||
      is.null(stats$xu) || is.null(stats$d))
    stop("Invalid 'stats' argument")
  if (length(stats$u) != n)
    stop("'stats$u' must be a vector of length 'n'")
  if (length(stats$xy) != p || length(stats$d) != p || length(stats$xu) != p)
    stop("'stats$xy', 'stats$xu' and 'stats$d' must be vectors of length 'p'")

  # Check inputs alpha0 and mu0.
  if (length(alpha0) != p || length(mu0) != p)
    stop("Inputs 'alpha0' and 'mu0' must be vectors of length 'p'")

  # Check input Xr0.
  if (length(Xr0) != n)
    stop("Input 'Xr0' must be a vector of length 'n'")

  # Check input S.
  if (sum(S < 1 | S > p) > 0)
      stop("Input 'S' contains invalid variable indices")

  # Execute the C routine, and return the results in a list object.
  # The only components of the list that change are alpha, mu and Xr.
  # We need to subtract 1 from the indices because R vector start at
  # 1, but C arrays start at 0. Note that I do not attempt to coerce X
  # here; if X is large, it could use a lot of memory to duplicate
  # this matrix. For this same reason, I set DUP = FALSE so that the
  # input arguments are not duplicated.
  result <- .C("varbvsbinupdateR",
               n       = as.integer(n),       # Number of samples.
               m       = as.integer(m),       # Number of updates.
               X       = X,                   # Matrix of samples.
               sa      = as.double(sa),       # Prior variance of coefficients.
               logodds = as.double(logodds),  # Prior log-odds.
               u       = as.double(stats$u),  # u = slope(eta).
               xy      = as.double(stats$xy), # xy = X'*yhat.
               xu      = as.double(stats$xu), # xu = X'*u.
               d       = as.double(stats$d),  # d = diag(X'*Uhat*X).
               alpha   = as.double(alpha0),   # Posterior inclusion prob's.
               mu      = as.double(mu0),      # Posterior mean coefficients.
               Xr      = as.double(Xr0),      # Xr = X*(alpha*mu).
               S       = as.integer(S-1),     # Updates to perform.
               DUP     = FALSE)
  return(list(alpha = result$alpha,
              mu    = result$mu,
              Xr    = result$Xr))
}
