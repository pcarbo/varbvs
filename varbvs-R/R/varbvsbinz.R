# This function implements the fully-factorized variational
# approximation for Bayesian variable selection in logistic
# regression, allowing for covariates. It is the same as varbvsbin,
# except that it allows for an additional set of covariates that are
# not subject to the same "spike-and-slab" priors as the other
# variables. The covariate data Z are specified as an n x m matrix,
# where n is the number of samples, and m is the number of
# covariates. This function is equivalent to varbvsbin when only one
# covariate is specified, the intercept, and Z = ones(n,1).
varbvsbinz <- function (X, Z, y, sa, logodds, alpha, mu, eta, tol = 1e-4,
                        maxiter = 1e4, verbose = TRUE, outer.iter = NULL,
                        update.sa = TRUE, optimize.eta = TRUE,n0 = 0,
                        sa0 = 0) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  # (1) INITIAL STEPS
  # -----------------
  # Compute a few useful quantities.
  Xr    <- c(X %*% (alpha*mu))
  stats <- updatestats_varbvsbinz(X,Z,y,eta)
  s     <- sa/(sa*stats$xdx + 1)

  # Initialize storage for outputs logw and err.
  logw <- rep(0,maxiter)
  err  <- rep(0,maxiter)
  
  # (2) MAIN LOOP
  # -------------
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  for (iter in 1:maxiter) {
    
    # Save the current variational parameters and model parameters.
    alpha0 <- alpha
    mu0    <- mu
    s0     <- s
    eta0   <- eta
    sa0    <- sa

    # (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    # --------------------------------------------
    logw0 <- int.logitz(Z,y,stats,alpha,mu,s,Xr,eta) +
             int.gamma(logodds,alpha) +
             int.klbeta(alpha,mu,s,sa)

    # (2b) UPDATE VARIATIONAL APPROXIMATION
    # -------------------------------------
    # Run a forward or backward pass of the coordinate ascent updates.
    if (iter %% 2)
      i <- 1:p
    else
      i <- p:1
    out   <- varbvsbinzupdate(X,sa,logodds,stats,alpha,mu,Xr,i)
    alpha <- out$alpha
    mu    <- out$mu
    Xr    <- out$Xr
    rm(out)

    # (2c) UPDATE ETA
    # ---------------
    # Update the free parameters specifying the variational approximation
    # to the logistic regression factors.
    if (optimize.eta) {
      eta   <- update_etaz(X,Z,y,betavar(alpha,mu,s),Xr,stats$d)
      stats <- updatestats_varbvsbinz(X,Z,y,eta)
      s     <- sa/(sa*stats$xdx + 1)
    }
    
    # (2d) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute variational lower bound to marginal log-likelihood.
    logw[iter] <- int.logitz(Z,y,stats,alpha,mu,s,Xr,eta) +
                  int.gamma(logodds,alpha) +
                  int.klbeta(alpha,mu,s,sa)

    # (2e) UPDATE PRIOR VARIANCE OF REGRESSION COEFFICIENTS
    # -----------------------------------------------------
    # Compute the maximum a posteriori estimate of sa, if requested.
    # Note that we must also recalculate the variance of the
    # regression coefficients when this parameter is updated.
    if (update.sa) {
      sa <- (sa0*n0 + dot(alpha,s + mu^2))/(n0 + sum(alpha))
      s  <- sa/(sa*stats$xdx + 1)
    }

    # (2f) CHECK CONVERGENCE
    # ----------------------
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum relative
    # difference between the parameters at two successive iterations
    # is less than the specified tolerance, or when the variational
    # lower bound has decreased. I ignore parameters that are very
    # small. If the variational bound decreases, stop.
    err[iter] <- max(abs(alpha - alpha0))
    if (verbose) {
      if (is.null(outer.iter))
        status <- NULL
      else
        status <- sprintf("%05d ",outer.iter)
      progress.str <- 
        paste(status,sprintf("%05d %+13.6e %0.1e %06.1f      NA %0.1e",
                             iter,logw[iter],err[iter],sum(alpha),sa),sep="")
      cat(progress.str)
      cat(rep("\r",nchar(progress.str)))
    }
    if (logw[iter] < logw0) {
      logw[iter]  <- logw0
      err[iter]   <- 0
      sa          <- sa0
      alpha       <- alpha0
      mu          <- mu0
      s           <- s0
      eta         <- eta0
      break
    } else if (err[iter] < tol)
      break
  }
    
  # Return the variational estimates.
  return(list(logw = logw[1:iter],err = err[1:iter],sa = sa,
              alpha = alpha,mu = mu,s = s,eta = eta))
}

# ----------------------------------------------------------------------
# diagprod(A,B) efficiently computes diag(A*B').
diagprod <- function (A, B)
  rowSums(A * B)

# ----------------------------------------------------------------------
# This function returns useful quantities for updating the variational
# approximation to the logistic regression factors, allowing for
# covariates.
updatestats_varbvsbinz <- function (X, Z, y, eta) {

  # Compute the slope of the conjugate.
  d <- slope(eta)
  
  # Compute the posterior covariance of u (coefficients for Z) given
  # beta (coefficients for X). This is equivalent to MATLAB expression
  # S = inv(Z'*diag(d)*Z).
  S <- solve(t(Z) %*% (Z * d))
  
  # Compute matrix dzr = D*Z*R', where D = diag(d), and R is an upper
  # triangular matrix such that R'*R = S.
  R   <- chol(S)
  dzr <- d * (Z %*% t(R))
  
  # Compute yhat. 
  yhat <- c(y - 0.5 - dzr %*% R %*% (t(Z) %*% (y - 0.5)))

  # Calculate xy = X'*yhat and xd = X'*d.
  xy <- c(yhat %*% X)
  xd <- c(d %*% X)
    
  # Compute the diagonal entries of X'*Dhat*X. For a definition of Dhat,
  # see the Bayesian Analysis journal paper.
  xdx <- diagsq(X,d) - diagsq(t(dzr) %*% X)

  # Return the result.
  return(list(S = S,d = d,yhat = yhat,xy = xy,xd = xd,xdx = xdx,dzr = dzr))
}

# ----------------------------------------------------------------------
# Computes the M-step update for the parameters specifying the
# variational lower bound to the logistic regression factors, allowing
# for additional covariates.
update_etaz <- function (X, Z, y, v, Xr, d) {

  # Compute the posterior covariance of u given beta. This is
  # equivalent to MATLAB expression S = inv(Z'*diag(d)*Z).
  S <- solve(t(Z) %*% (Z * d))

  # Compute the posterior mean of the regression coefficients
  # corresponding to the covariates.
  muz <- S %*% t(Z) %*% (y - 0.5 - d*Xr)
  
  # Calculate the covariance between the coefficients u and beta.
  W <- (-t((S %*% t(Z * d)) %*% X) * v)
  
  # This is the M-step update for the free parameters.
  U <- t((S %*% t(Z * d)) %*% X)
  return(sqrt(c(Z %*% muz + Xr)^2 + diagsq2(Z,S) +
              diagsq2(Z,t(U) %*% (U * v)) + diagsqt(X,v) +
              2*diagprod(X %*% W,Z)))
}

# -----------------------------------------------------------------------
int.logitz <- function (Z, y, stats, alpha, mu, s, Xr, eta) {

  # Get some of the statistics.
  yhat <- stats$yhat
  xdx  <- stats$xdx
  S    <- stats$S
  d    <- stats$d

  # Compute the variational approximation to the expectation of the
  # log-likelihood with respect to the approximate posterior distribution.
  return(sum(logsigmoid(eta)) + dot(eta,d*eta - 1)/2 +
         c(determinant(S,logarithm = TRUE)$modulus/2) +
         qnorm(t(Z) %*% (y - 0.5),S)^2/2 + dot(yhat,Xr) - qnorm(Xr,d)^2/2 +
         qnorm(t(Z) %*% (Xr * d),S)^2/2 - dot(xdx,betavar(alpha,mu,s))/2)
}
