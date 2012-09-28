varbvsoptimize <- function (X, y, sigma, sa, logodds, alpha0 = NULL,
                            mu0 = NULL, verbose = TRUE) {
  # Implements the fully-factorized variational approximation for
  # Bayesian variable selection in linear regression. It finds the
  # "best" fully-factorized variational approximation to the posterior
  # distribution of the coefficients in a linear regression model of a
  # continuous outcome (quantitiative trait), with spike and slab
  # priors on the coefficients. By "best", we mean the approximating
  # distribution that locally minimizes the Kullback-Leibler
  # divergence between the approximating distribution and the exact
  # posterior.

  # Convergence is reached when the maximum relative distance between
  # successive updates of the variational parameters is less than this
  # quantity.
  tolerance <- 1e-4  

  # CHECK INPUTS.
  # Check input X.
  if (!is.matrix(X))
    stop("Input argument 'X' must be a matrix")
  if (!is.double(X))
    storage.mode(X) <- "double"

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  # Check input y.
  y <- c(y)
  if (length(y) != n)
    stop("Data 'X' and 'y' do not match")

  # Check inputs sigma and sa.
  if (!is.scalar(sigma) || !is.scalar(sa))
    stop("Input arguments 'sigma' and 'sa' must be scalars")

  # Check input logodds.
  if (is.scalar(logodds))
    logodds <- rep(logodds,p)
  if (length(logodds) != p)
    stop("Input 'logodds' must be a scalar or a vector of length 'p'")

  # TAKE CARE OF OPTIONAL INPUTS.
  # Set initial estimates of variational parameters.
  if (is.null(alpha0)) {
    alpha <- runif(p)
    alpha <- alpha / sum(alpha)
  }
  else
    alpha <- c(alpha0)
  if (is.null(mu0))
    mu <- rnorm(p)
  else
    mu <- c(mu0)
  if (length(alpha) != p || length(mu) != p)
    stop("'alpha0' and 'mu0' must be NULL or vectors of length 'p'")
  
  # INITIAL STEPS.
  # Compute a few useful quantities.
  xy <- c(y %*% X)           # xy = X'*y.
  d  <- diagsq(X)            # d  = diag(X'*Y).
  Xr <- c(X %*% (alpha*mu))  # Xr = X*(alpha*mu).

  # Calculate the variance of the coefficients.
  s <- sa*sigma/(sa*d + 1)

  # MAIN LOOP.
  # Repeat until convergence criterion is met.
  lnZ  <- -Inf
  iter <- 0
  if (verbose) {
    cat("       variational    max. incl max.\n")
    cat("iter   lower bound  change vars E[b]\n")
  }
  while (TRUE) {

    # Go to the next iteration.
    iter <- iter + 1

    # Save the current variational parameters and lower bound.
    alpha0  <- alpha
    mu0     <- mu
    lnZ0    <- lnZ
    params0 <- c(alpha,alpha*mu)

    # UPDATE VARIATIONAL APPROXIMATION.
    # Run a forward or backward pass of the coordinate ascent updates.
    if (is.odd(iter))
      S <- seq(1,p)
    else
      S <- seq(p,1,-1)
    result <- varbvsupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,S)
    alpha  <- result$alpha
    mu     <- result$mu
    Xr     <- result$Xr

    # COMPUTE VARIATIONAL LOWER BOUND.
    # Compute the lower bound to the marginal log-likelihood.
    lnZ <- intlinear(Xr,d,y,sigma,alpha,mu,s) +
           intgamma(logodds,alpha) +
           intklbeta(alpha,mu,s,sigma*sa)

    # CHECK CONVERGENCE.
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum relative
    # difference between the parameters at two successive iterations
    # is less than the specified tolerance, or when the variational
    # lower bound has decreased. I ignore parameters that are very
    # small.
    params <- c(alpha,alpha*mu)
    S      <- which(abs(params) > 1e-6)
    err    <- relerr(params[S],params0[S])
    if (verbose)
      cat(sprintf('%4d %+13.6e %0.1e %4d %0.2f\n',iter,lnZ,max(err),
                  round(sum(alpha)),max(abs(alpha*mu))))
    if (lnZ < lnZ0) {
      alpha <- alpha0
      mu    <- mu0
      lnZ   <- lnZ0
      break
    }
    else if (max(err) < tolerance)
      break
  }

  # Return the variational estimates.
  return(list(alpha=alpha,mu=mu,s=s,lnZ=lnZ))
}
