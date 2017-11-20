# TO DO: Explain what this function does, and how to use it.
varbvssparse <- function (X, y, k, sigma, sa, logodds, alpha, mu, tol = 1e-4,
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
    cat("        variational    max.\n")
    cat(" iter   lower bound  change\n")
  }
  for (iter in 1:maxiter) {

    # Save the current variational parameters.
    alpha0 <- alpha
    mu0    <- mu
    
    # (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood (the
    # "ELBO").
    logw0 <- varbvssparse.elbo(y,sigma,sa,alpha,mu,Xr)

    # (2b) UPDATE VARIATIONAL APPROXIMATION
    # -------------------------------------
    # Run a forward or backward pass of the coordinate ascent updates.
    if (iter %% 2)
      i <- 1:k
    else
      i <- k:1
    out <- varbvssparseupdate(X,sigma,sa,logodds,xy,s,alpha,mu,Xr,i)
    alpha <- out$alpha
    mu    <- out$mu
    Xr    <- out$Xr
    rm(out)

    # (2c) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood (the
    # "ELBO").
    logw[iter] <- varbvssparse.elbo(y,sigma,sa,alpha,mu,Xr)

    # (2d) CHECK CONVERGENCE
    # ----------------------
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum difference
    # between the posterior probabilities at two successive iterations
    # is less than the specified tolerance, or when the variational
    # lower bound has decreased.
    err[iter] <- max(abs(alpha - alpha0))
    if (verbose) {
      progress.str <- sprintf("%05d %+13.6e %0.1e",iter,logw[iter],err[iter])
      cat(progress.str)
      cat(rep("\r",nchar(progress.str)))
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

# Execute a single round of the coordinate ascent updates to maximize
# the variational lower bound for the "sparse" Bayesian variable
# selection model.
varbvssparseupdate <- function (X, sigma, sa, logodds, xy, s,
                                alpha0, mu0, Xr0, i) {

  # Initialize storage for the results.
  alpha <- alpha0
  mu    <- mu0
  Xr    <- Xr0

  # Repeat for each effect to update
  for (j in i) {
    
    # Update the variational estimate of the posterior means.
    mu[,j] <- s[,j]/sigma * (xy - c(rowSums(Xr[,-j]) %*% X))
    
    # Update the variational estimate of the posterior inclusion
    # probabilities.
    alpha[,j] <- sigmoid(logodds + (log(s[,j]/(sa*sigma)) + mu[,j]^2/s[,j])/2)
    alpha[,j] <- alpha[,j]/sum(alpha[,j])
    
    # Update the ith effect.
    Xr[,j] <- X %*% (alpha[,j]*mu[,j])
  }
  
  return(list(alpha = alpha,mu = mu,Xr = Xr))
}

# TO DO: Explain here what this function does, and how to use it.
varbvssparse.elbo <- function (y, sigma, sa, alpha, mu, Xr) {
  return(- length(y)/2*log(2*pi*sigma)
         - norm2(y - rowSums(Xr))^2/(2*sigma)
         + sum(apply(Xr,2,norm2)^2)/(2*sigma)
         - sum(d %*% betavar(alpha,mu,s))/(2*sigma)
         + sum(alpha*logsigmoid(logodds))
         - sum(alpha*log(alpha + eps))
         + sum(alpha*(1 + log(s/sa) - (mu^2 + s)/sa)/2))
}

# Replicate vector x to create an n x m matrix, where m = length(x).
rep.row <- function (x, n)
  matrix(x,n,length(x),byrow = TRUE)

# Retrieve a couple hidden functions from the varbvs package.
sigmoid    <- varbvs:::sigmoid
logsigmoid <- varbvs:::logsigmoid
diagsq     <- varbvs:::diagsq
norm2      <- varbvs:::norm2
  
# Shorthand for machine precision.
eps <- .Machine$double.eps
