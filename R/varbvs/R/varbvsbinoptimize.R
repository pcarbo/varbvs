varbvsbinoptimize <- function (X, y, sa, logodds, alpha0 = NULL, mu0 = NULL,
                               eta0 = NULL, fixed.eta = FALSE,
                               verbose = TRUE) {
  # Implements the fully-factorized variational approximation for
  # Bayesian variable selection in logistic regression. It finds the
  # "best" fully-factorized variational approximation to the posterior
  # distribution of the coefficients in a logistic regression model of
  # a binary outcome or trait (e.g. disease status in a case-control
  # study), with spike and slab priors on the coefficients. By "best",
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
  if (!is.scalar(sa))
    stop("Input argument 'sa' must be a scalar")

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

  # Determine whether to update the variational approximation to the
  # logistic regression.
  if (fixed.eta && is.null(eta0))
    stop("'fixed.eta = TRUE' requires specification of argument 'eta0'")
  
  # Initialize the free parameters specifying the variational
  # approximation to the logistic regression factors.
  if (is.null(eta0))
    eta <- rep(1,n)
  else
    eta <- c(eta0)
  if (length(eta) != n)
    stop("Input argument 'eta' must NULL or a vector of length 'n'")

  # INITIAL STEPS.
  # Compute a few useful quantities.
  Xr    <- c(X %*% (alpha*mu))  # Xr = X*(alpha*mu).
  stats <- updatestats(X,y,eta)

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
    eta0    <- eta
    params0 <- c(alpha,alpha*mu)

    # UPDATE VARIATIONAL APPROXIMATION.
    # Run a forward or backward pass of the coordinate ascent updates.
    if (is.odd(iter))
      S <- seq(1,p)
    else
      S <- seq(p,1,-1)
    result <- varbvsbinupdate(X,sa,logodds,stats,alpha,mu,Xr,S)
    alpha  <- result$alpha
    mu     <- result$mu
    Xr     <- result$Xr

    # Recalculate the posterior variance of the coefficients.
    s <- sa/(sa*stats$d + 1)

    # UPDATE ETA.
    # Update the free parameters specifying the variational approximation
    # to the logistic regression factors.
    if (!fixed.eta) {
      eta   <- update.eta(X,y,betavar(alpha,mu,s),Xr,stats$u);
      stats <- updatestats(X,y,eta)
    }

    # COMPUTE VARIATIONAL LOWER BOUND.
    # Compute variational lower bound to marginal log-likelihood.
    lnZ <- intlogit(y,stats,alpha,mu,s,Xr,eta) +
           intgamma(logodds,alpha) +
           intklbeta(alpha,mu,s,sa)

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
      eta   <- eta0
      lnZ   <- lnZ0
      break
    }
    else if (max(err) < tolerance)
      break
  }

  # Return the variational estimates.
  return(list(alpha=alpha,mu=mu,s=s,eta=eta,lnZ=lnZ))
}

slope <- function (x) {
  # Compute (sigmoid(x) - 1/2)/x, the slope of the conjugate to the
  # log-sigmoid function at x, times 2. For details, see Bishop
  # (2006), or the Bayesian Analysis paper. This is useful for working
  # with the variational approximation to the logistic regression.
  return((sigmoid(x) - 0.5)/(x + eps))
}

updatestats <- function (X, y, eta) {
  # Calculate some useful quantities for updating the variational
  # approximation to the logistic regression factors. Inputs X and y
  # specify the data, and eta is the vector of free parameters. It is
  # a column vector of length equal to the number of samples. This
  # function should be called whenever the free parameters are
  # modified.

  # Compute the slope of the conjugate.
  u <- slope(eta)

  # Compute 'beta0' and 'yhat'. See the journal paper for an
  # explanation of these two variables.
  beta0 <- sum(y - 0.5)/sum(u)
  yhat  <- y - 0.5 - beta0*u

  # Calculate a couple useful matrix-vector products.
  xy <- c(y %*% X)  # xy = X'*y.
  xu <- c(u %*% X)  # xu = X'*u.

  # Compute the diagonal entries of X'*Uhat*X. For a definition of
  # matrix Uhat, see the journal paper.
  d <- diagsq(X,u) - xu^2/sum(u)

  # Return the quantities as a list.
  return(list(u=u,yhat=yhat,xy=xy,xu=xu,d=d))
}

update.eta <- function (X, y, v, Xr, u) {
  # Returns the M-step update for the parameters specifying the
  # variational lower bound to the logistic regression factors. Input
  # Xr be equal to Xr = X*(alpha*mu). Input v is the posterior
  # variance of the coefficients. For this update to be valid, the
  # posterior covariance of the coefficients must equal to diag(v).
  # Input u must be u = slope(eta), where u is the current estimate of
  # the free parameters; see function 'slope' for details.
  # 
  # Note that vector y and matrix X must *not* be centered.  
  
  # Compute 'mu0', the posterior mean of the intercept in the logistic
  # regression under the variational approximation. Here, 'a' is the
  # variance of the intercept given the other coefficients.
  a   <- 1/sum(u)
  mu0 <- a*(sum(y - 0.5) - dot(u,Xr))

  # Compute 's0', the (marginal) posterior variance of the intercept in
  # the logistic regression.
  xu <- c(u %*% X)  # xu = X'*u.
  s0 <- a*(1 + a*dot(v,xu^2))
  
  # Calculate the covariance between the intercept and coefficients.
  c0 <- -a*xu*v
  
  # This is the M-step update for the free parameters.
  eta <- sqrt((mu0 + Xr)^2 + s0 + diagsqt(X,v) + 2*c(X %*% c0))
  return(eta)
}

intlogit <- function (y, stats, alpha, mu, s, Xr, eta) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood for the logistic regression
  # model. This integral is an approximation to the expectation of the
  # logistic regression log-likelihood taken with respect to the
  # variational approximation. Input argument 'stats' must be the
  # return argument to updatestats(X,y,eta). Input Xr be equal to Xr
  # = X*(alpha*mu).

  # Get the variance of the intercept given the other coefficients.
  u <- stats$u
  a <- 1/sum(u)

  # Compute the variational approximation to the expectation of the
  # log-likelihood with respect to the variational approximation.
  f <- sum(logsigmoid(eta)) + dot(eta,u*eta - 1)/2 + log(a)/2 +
       a*sum(y - 0.5)^2/2 + dot(stats$yhat,Xr) - qnorm(Xr,u)^2/2 +
       a/2*dot(u,Xr)^2 - dot(stats$d,betavar(alpha,mu,s))/2
  return(f)
}
