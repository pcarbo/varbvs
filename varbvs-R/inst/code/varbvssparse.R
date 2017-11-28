# Shorthand for machine precision.
eps <- .Machine$double.eps

# Retrieve a couple hidden functions from the varbvs package.
rep.row    <- varbvs:::rep.row
rep.col    <- varbvs:::rep.col
sigmoid    <- varbvs:::sigmoid
logsigmoid <- varbvs:::logsigmoid
diagsq     <- varbvs:::diagsq
norm2      <- varbvs:::norm2
betavar    <- varbvs:::betavar

# Implements the Sum of Single Effects (SuSiE) model and variational
# approximation. Input arguments include:
#
#   X is an n x p matrix of observations of the variables (or
#   features), where n is the number of samples, and p is the number
#   of variables
#
#   y is a vector containing samples of the continuous outcome; it is
#   a vector of length n
#
#   k > 1 specifies the number of (nonzero) effects in the linear
#   regression (the "trivial" case of k = 1 is not yet implemented)
#
#   sigma is a scalar giving the variance of the residual
#
#   sa is a scalar specifying the prior variance of the regression
#   coefficients (this prior variance is scaled by sigma)
#
#   pp specifies the prior probability that a coefficient is nonzero;
#   it is a vector of length p
#
#   alpha is an optional argument giving initial estimates of the
#   variational parameter "alpha"; it should be an p x k matrix
#
#   mu is an optional argument giving initial estimates of the
#   variational parameter "mu"; it should be a p x k matrix
#
varbvssparse <- function (X, y, k, sigma, sa, pp, alpha, mu, tol = 1e-4,
                          maxiter = 1e4, verbose = TRUE) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  # Check input "k".
  if (k < 2)
    stop("Argument \"k\" should be at least 2.")
  
  # Set initial estimates of variational parameter alpha.
  if (missing(alpha)) {
    alpha <- rand(p,k)
    alpha <- alpha / rep.row(colSums(alpha),p)
  }

  # Set initial estimates of variational parameter mu.
  if (missing(mu))
    mu <- randn(p,k)
  
  # (1) INITIAL STEPS
  # -----------------
  # Compute a few useful quantities. 
  xy <- c(y %*% X)
  d  <- diagsq(X)
  Xr <- X %*% (alpha*mu)

  # Calculate the variance of the coefficients.
  s <- sa*sigma/(sa*d + 1)
  s <- matrix(s,p,k)
  
  # Initialize storage for outputs logw and err.
  logw <- rep(0,maxiter)
  err  <- rep(0,maxiter)
  
  # (2) MAIN LOOP
  # -------------
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  if (verbose) {
    cat("             variational    max.\n")
    cat("iter         lower bound  change\n")
  }
  for (iter in 1:maxiter) {

    # Save the current variational parameters.
    alpha0 <- alpha
    mu0    <- mu
    
    # (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood ("ELBO").
    logw0 <- varbvssparse.elbo(y,d,sigma,sigma*sa,pp,alpha,mu,s,Xr)

    # (2b) UPDATE VARIATIONAL APPROXIMATION
    # -------------------------------------
    # Run a forward or backward pass of the coordinate ascent updates.
    if (iter %% 2)
      i <- 1:k
    else
      i <- k:1
    out   <- varbvssparseupdate(X,sigma,sa,pp,xy,s,alpha,mu,Xr,i)
    alpha <- out$alpha
    mu    <- out$mu
    Xr    <- out$Xr
    rm(out)
    
    # (2c) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood ("ELBO").
    logw[iter] <- varbvssparse.elbo(y,d,sigma,sigma*sa,pp,alpha,mu,s,Xr)

    # (2d) CHECK CONVERGENCE
    # ----------------------
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum difference
    # between the posterior probabilities at two successive iterations
    # is less than the specified tolerance, or when the variational
    # lower bound has decreased.
    err[iter] <- max(abs(alpha - alpha0))
    if (verbose) {
      progress.str <- sprintf("%04d %+13.12e %0.1e",iter,logw[iter],err[iter])
      cat(progress.str)
      cat("\n")
    }
    if (logw[iter] < logw0) {
      err[iter]  <- 0
      logw[iter] <- logw0
      alpha      <- alpha0
      mu         <- mu0
      break
    } else if (err[iter] < tol)
      break
  }
  if (verbose)
    cat("\n")
  return(list(logw = logw[1:iter],err = err[1:iter],sigma = sigma,sa = sa,
              alpha = alpha,mu = mu,s = s))
}

# ----------------------------------------------------------------------
# Execute a single round of the coordinate ascent updates to maximize
# the variational lower bound for the "sparse" Bayesian variable
# selection model.
varbvssparseupdate <- function (X, sigma, sa, pp, xy, s,
                                alpha0, mu0, Xr0, i) {

  # Initialize storage for the results.
  alpha <- alpha0
  mu    <- mu0
  Xr    <- Xr0

  # Repeat for each effect to update.
  for (j in i) {
    
    # Update the variational estimate of the posterior means.
    mu[,j] <- s[,j]/sigma * (xy - c(rowSums(Xr[,-j]) %*% X))
    
    # Update the variational estimate of the posterior inclusion
    # probabilities.
    alpha[,j] <- log(pp + eps) + (log(s[,j]/(sa*sigma)) + mu[,j]^2/s[,j])/2
    alpha[,j] <- normalizelogweights(alpha[,j])
    
    # Update the jth set of effects.
    Xr[,j] <- X %*% (alpha[,j]*mu[,j])
  }
  
  return(list(alpha = alpha,mu = mu,Xr = Xr))
}

# ----------------------------------------------------------------------
# Compute the variational lower bound to the marginal log-likelihood
# ("ELBO").
varbvssparse.elbo <- function (y, d, sigma, sa, pp, alpha, mu, s, Xr)
  (-length(y)/2*log(2*pi*sigma)
   - norm2(y - rowSums(Xr))^2/(2*sigma) + sum(Xr^2)/(2*sigma)
   - sum(d %*% (alpha*(mu^2 + s)))/(2*sigma)
   + sum(log(pp + eps) %*% alpha) - sum(alpha*log(alpha + eps))
   + sum(alpha*(1 + log(s/sa) - (mu^2 + s)/sa)/2))

