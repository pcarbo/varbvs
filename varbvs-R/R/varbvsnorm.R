# Implements the fully-factorized variational approximation for
# Bayesian variable selection in linear regression. It finds the
# "best" fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of a
# continuous outcome (quantitiative trait), with spike and slab priors
# on the coefficients. By "best", we mean the approximating
# distribution that locally minimizes the K-L divergence between the
# approximating distribution and the exact posterior.
#
# Input X is an n x p matrix of observations about the variables (or
# features), where n is the number of samples, and p is the number of
# variables. Input y is the vector of observations about the outcome;
# it is a vector of length n.
#
# Inputs sigma, sa and logodds are additional model parameters; sigma
# and sa are scalars. Input sigma specifies the variance of the
# residual, and sa specifies the prior variance of the
# coefficients. logodds is the prior log-odds of inclusion for each
# variable.
#
# Output logw is the variational estimate of the marginal
# log-likelihood given the hyperparameters. Outputs alpha, mu and s
# are the parameters of the variational approximation and,
# equivalently, variational estimates of posterior quantites: under
# the variational approximation, the ith regression coefficient is
# normal with probability alpha(i); mu(i) and s(i) are the mean and
# variance of the coefficient given that it is included in the model.
#
# When update.sa = TRUE, there is the additional option of computing
# the maximum a posteriori (MAP) estimate of the prior variance
# parameter (sa), in which sa is drawn from a scaled inverse
# chi-square distribution with scale sa0 and degrees of freedom n0.
varbvsnorm <-
  function (X, y, sigma, sa, logodds, alpha, mu, tol = 1e-4, maxiter = 1e4,
            verbose = TRUE, outer.iter = NULL, update.sigma = TRUE,
            update.sa = TRUE, n0 = 0, sa0 = 0) {

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)
  
  # (1) INITIAL STEPS
  # -----------------
  # Compute a few useful quantities. 
  xy <- c(y %*% X)
  d  <- diagsq(X)
  Xr <- c(X %*% (alpha*mu))
  
  # Calculate the variance of the coefficients.
  s <- sa*sigma/(sa*d + 1)

  # (2) MAIN LOOP
  # -------------
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  logw <- (-Inf)
  for (iter in 1:maxiter) {

    # Save the current variational and model parameters.
    alpha0 <- alpha
    mu0    <- mu
    sigma0 <- sigma

    # (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
    logw0 <- int.linear(Xr,d,y,sigma,alpha,mu,s) +
             int.gamma(logodds,alpha) +
             int.klbeta(alpha,mu,s,sigma*sa)

    # (2b) UPDATE VARIATIONAL APPROXIMATION
    # -------------------------------------
    # Run a forward or backward pass of the coordinate ascent updates.
    if (iter %% 2)
      i <- 1:p
    else
      i <- p:1
    out   <- varbvsnormupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,i)
    alpha <- out$alpha
    mu    <- out$mu
    Xr    <- out$Xr
    rm(out)

    browser()
    
    # (2c) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
    logw <- int.linear(Xr,d,y,sigma,alpha,mu,s) +
            int.gamma(logodds,alpha) +
            int.klbeta(alpha,mu,s,sigma*sa)

    # (2d) UPDATE RESIDUAL VARIANCE
    # -----------------------------
    # Compute the maximum likelihood estimate of sigma, if requested.
    # Note that we must also recalculate the variance of the regression
    # coefficients when this parameter is updated.
    if (update.sigma) {
      sigma <- (norm2(y - Xr)^2 + dot(d,betavar(alpha,mu,s)) +
                dot(alpha,(s + mu^2)/sa))/(n + sum(alpha))
      s     <- sa*sigma/(sa*d + 1)
    }
    
    # (2e) UPDATE PRIOR VARIANCE OF REGRESSION COEFFICIENTS
    # -----------------------------------------------------
    # Compute the maximum a posteriori estimate of sa, if requested.
    # Note that we must also recalculate the variance of the
    # regression coefficients when this parameter is updated.
    if (update.sa) {
      sa <- (sa0*n0 + dot(alpha,s + mu^2))/(n0 + sigma*sum(alpha))
      s  <- sa*sigma/(sa*d + 1)
    }

    # (2f) CHECK CONVERGENCE
    # ----------------------
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum difference
    # between the posterior probabilities at two successive iterations
    # is less than the specified tolerance, or when the variational
    # lower bound has decreased.
    err <- abs(alpha - alpha0)
    if (verbose) {
      if (isempty(outer.iter))
        status <- ""
      else
        status <- sprintf("%05d ",outer.iter)
      status <- paste(status,sprintf("%05d %+13.6e %0.1e %06.1f %0.1e %0.1e",
                                     iter,logw,max(err),sum(alpha),
                                     sigma,sa),sep="")
      caterase(status)
    }
    if (logw < logw0) {
      return(list(logw = logw0,sigma = sigma0,sa = sa0,
                  alpha = alpha0,mu = mu0,s = s))
      logw  <- logw0
      sa    <- sa0
      alpha <- alpha0
      mu    <- mu0
      break
    } else if (max(err) < tol)
      return(list(logw = logw,sigma = sigma,sa = sa,
                  alpha = alpha,mu = mu,s = s))
  }
}

# ----------------------------------------------------------------------
# Computes an integral that appears in the variational lower bound of
# the marginal log-likelihood. This integral is the expectation of the
# linear regression log-likelihood taken with respect to the
# variational approximation.
int.linear <- function (Xr, d, y, sigma, alpha, mu, s) {
  n <- length(y)
  return(-length(y)/2*log(2*pi*sigma) - norm2(y - Xr)^2/(2*sigma) 
         - dot(d,betavar(alpha,mu,s))/(2*sigma))
}
