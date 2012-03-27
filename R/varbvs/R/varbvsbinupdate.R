## % [ALPHA,MU,XR] = VARBVSBINUPDATE(X,SA,LOGODDS,STATS,ALPHA0,MU0,XR0,I) runs
## % a single iteration of the coordinate ascent updates to maximize the
## % variational lower bound for Bayesian variable selection in logistic
## % regression. It adjusts the fully-factorized variational approximation to
## % the posterior distribution of the coefficients in a logistic regression
## % model of a binary outcome or trait, with spike and slab priors on the
## % coefficients.
## %
## % All inputs to this function are required. Input X is an N x P matrix of
## % observations about the variables (or features), where N is the number of
## % samples, and P is the number of variables. Y is the vector of observations
## % about the binary trait; it is a vector of length N. Unlike function
## % VARBVSUPDATE, Y and X must *not* be centered. Instead, we will account for
## % the intercept as we update the variational approximation.
## %
## % This routine is implemented with the assumption that X is a single
## % floating-point precision matrix (type HELP SINGLE), as opposed to MATLAB's
## % default of double precision. This is useful for large data sets, because
## % single precision requires half of the number of bits as double
## % floating-point precision. If X is provided in another numerical
## % representation, an error is reported.
## %
## % Input scalar SA specifies the prior variance of the coefficients. LOGODDS
## % is the prior log-odds of inclusion for each variable. It is equal to
## % LOGODDS = LOG(Q./(1-Q)), where Q is the prior probability that each
## % variable is included in the linear model of Y. LOGODDS is a vector of
## % length P. Note that a residual variance parameter SIGMA is not needed to
## % model a binary trait. Input STATS is the STRUCT output from UPDATESTATS.
## %
## % Inputs ALPHA0, MU0 are the current parameters of the variational
## % approximation; under the variational approximation, the ith regression
## % coefficient is normal with probability ALPHA0(i), and MU0(i) is the mean
## % of the coefficient given that it is included in the model. Inputs XR0 must
## % be XR0 = X*(ALPHA0.*MU0).
## %
## % Input I specifies the order in which the coordinates are updated. It may
## % be a vector of any length. Each entry of I must be an integer between 1
## % and P.
## %
## % There are three outputs. Output vectors ALPHA and MU are the updated
## % variational parameters, and XR = X*(ALPHA.*MU). The computational
## % complexity of VARBVSBINUPDATE is O(N*LENGTH(I)).
varbvsbinupdate <- function (X, sa, logodds, stats, alpha0, mu0, Xr0, S) {
  # *** TO DO: BRIEF DESCRIPTION OF FUNCTION GOES HERE. ***

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
