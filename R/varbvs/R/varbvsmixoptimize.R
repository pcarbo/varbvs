varbvsmixoptimize <- function (X, y, sigma, sa1, sa2, logodds, alpha0 = NULL,
                               mu10 = NULL, mu20 = NULL, verbose = TRUE) {
  # Implements the fully-factorized variational approximation for the
  # Bayesian variable selection mixture model with linear
  # regression. It finds the "best" fully-factorized variational
  # approximation to the posterior distribution of the coefficients in
  # a linear regression model of a continuous outcome (quantitiative
  # trait), with mixture of normal priors on the coefficients. By "best",
  # we mean the approximating distribution that locally minimizes the
  # Kullback-Leibler divergence between the approximating distribution
  # and the exact posterior.

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
  if (!is.scalar(sigma) || !is.scalar(sa1) || !is.scalar(sa2))
    stop("Input arguments 'sigma', 'sa1' and 'sa2' must be scalars")

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
  if (is.null(mu10))
    mu1 <- sqrt(sigma*sa1) * rnorm(p)
  else
    mu1 <- c(mu10)
  if (is.null(mu20))
    mu2 <- sqrt(sigma*sa2) * rnorm(p)
  else
    mu2 <- c(mu20)
  if (length(alpha) != p || length(mu1) != p || length(mu2) != p)
    stop("'alpha0', 'mu10' and 'mu20' must be NULL or vectors of length 'p'")
  
  # INITIAL STEPS.
  # Compute a few useful quantities.
  xy <- c(y %*% X)  # xy = X'*y.
  d  <- diagsq(X)   # d  = diag(X'*X).
  Xr <- c(X %*% betameanmix(alpha,mu1,mu2))

  # Calculate the variance of the coefficients.
  s1 <- sa1*sigma/(sa1*d + 1)
  s2 <- sa2*sigma/(sa2*d + 1)

  # MAIN LOOP.
  # Repeat until convergence criterion is met.
  lnZ  <- -Inf
  iter <- 0
  if (verbose) {
    cat("       variational    max. incl max.      \n")
    cat("iter   lower bound  change vars E[b] sigma\n")
  }
  while (TRUE) {

    # Go to the next iteration.
    iter <- iter + 1

    # Save the current variational parameters and lower bound.
    alpha0  <- alpha
    mu10    <- mu1
    mu20    <- mu2
    lnZ0    <- lnZ
    params0 <- c(alpha,alpha*mu1,(1-alpha)*mu2)

    # UPDATE VARIATIONAL APPROXIMATION.
    # Run a forward or backward pass of the coordinate ascent updates.
    if (is.odd(iter))
      S <- seq(1,p)
    else
      S <- seq(p,1,-1)
    out   <- varbvsmixupdate(X,sigma,sa1,sa2,logodds,xy,d,alpha,mu1,mu2,Xr,S)
    alpha <- out$alpha
    mu1   <- out$mu1
    mu2   <- out$mu2
    Xr    <- out$Xr

    # UPDATE RESIDUAL VARIANCE.
    # Compute the maximum likelihood estimate of the residual variance
    # parameter, sigma. Note that after updating the residual variance
    # parameter, we must also recalculate the variance of the regression
    # coefficients. 
    sigma <- (norm2(y - Xr)^2 + dot(d,betavarmix(alpha,mu1,mu2,s1,s2))
              + dot(alpha,s1 + mu1^2)/sa1
              + dot(1-alpha,s2 + mu2^2)/sa2)/(n + p)
    s1    <- sa1*sigma/(sa1*d + 1)
    s2    <- sa2*sigma/(sa2*d + 1)
    
    # COMPUTE VARIATIONAL LOWER BOUND.
    # Compute the lower bound to the marginal log-likelihood.
    lnZ <- (intlinearmix(Xr,d,y,sigma,alpha,mu1,mu2,s1,s2) 
            + intgamma(logodds,alpha) 
            + intklbetamix(alpha,mu1,mu2,s1,s2,sigma*sa1,sigma*sa2))

    # CHECK CONVERGENCE.
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum relative
    # difference between the parameters at two successive iterations
    # is less than the specified tolerance, or when the variational
    # lower bound has decreased. I ignore parameters that are very
    # small.
    params <- c(alpha,alpha*mu1,(1-alpha)*mu2)
    S      <- which(abs(params) > 1e-6)
    err    <- relerr(params[S],params0[S])
    if (verbose)
      cat(sprintf('%4d %+13.6e %0.1e %4d %0.2f %5.2f\n',iter,lnZ,max(err),
                  round(sum(alpha)),max(abs(alpha*mu1)),sqrt(sigma)))
    if (lnZ < lnZ0) {
      alpha <- alpha0
      mu1   <- mu10
      mu2   <- mu20
      lnZ   <- lnZ0
      break
    }
    if (max(err) < tolerance)
      break
  }

  # Return the variational estimates.
  return(list(lnZ=lnZ,alpha=alpha,mu1=mu1,mu2=mu2,s1=s1,s2=s2,sigma=sigma))
}
