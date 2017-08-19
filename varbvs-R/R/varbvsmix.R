# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2017, Peter Carbonetto
#
# This program is free software: you can redistribute it under the
# terms of the GNU General Public License; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# Fit linear regression model with mixture-of-normals prior using 
# variational approximation techniques. See varbvsmix.Rd for details.
varbvsmix <- function (X, Z, y, sa, sigma, w, alpha, mu, update.sigma,
                       update.sa, update.w, w.penalty, drop.threshold = 1e-8,
                       tol = 1e-4, maxiter = 1e4, verbose = TRUE) {
    
  # Get the number of samples (n), variables (p) and mixture
  # components (K).
  n <- nrow(X)
  p <- ncol(X)
  K <- length(sa)

  # (1) CHECK INPUTS
  # ----------------
  # Check input X.
  if (!(is.matrix(X) & is.double(X) & sum(is.na(X)) == 0))
    stop("Input X must be a double-precision matrix with no missing values.")

  # Add column names to X if they are not provided.
  if (is.null(colnames(X)))
    colnames(X) <- paste0("X",1:p)

  # Check input Z.
  if (!is.null(Z)) {
    Z <- as.matrix(Z)
    if (!is.numeric(Z) | sum(is.na(Z)) > 0)
      stop("Input Z must be a numeric matrix with no missing values.")
    if (nrow(Z) != n)
      stop("Inputs X and Z do not match.")
    storage.mode(Z) <- "double"
  }

  # Add intercept.
  if (is.null(Z))
    Z <- matrix(1,n,1)
  else
    Z <- cbind(1,Z)

  # Check input y.
  if (!is.numeric(y) | sum(is.na(y)) > 0)
    stop("Input y must be a numeric vector with no missing values.")
  if (length(y) != n)
    stop("Inputs X and y do not match")
  y <- c(as.double(y))
  
  # Check input sa. The variance of the first mixture component should
  # be exactly zero, corresponding to the "spike".
  sa <- c(as.double(sa))
  if (sa[1] != 0)
    stop("Variance of first mixture component should be 0")

  # (2) PROCESS OPTIONS
  # -------------------
  if (!is.finite(maxiter))
    stop("Input maxiter must be a finite number")
  
  # Get initial estimate for the variance of the residual, if
  # provided.
  if (missing(sigma)) {
    sigma <- var(y)
    update.sigma.default <- TRUE
  } else {
    sigma <- c(sigma)
    update.sigma.default <- FALSE
  }

  # Get initial estimate for the mixture weights, if provided.
  if (missing(w))
    w <- rep(1/K,K)
  else
    w <- c(w)
  if (length(w) != K)
    stop("Length of input w should be the same as length of sa")
  
  # Specify the penalty term for the mixture weights.
  if (missing(w.penalty))
    w.penalty <- rep(1,K)
  else
    w.penalty <- c(w.penalty)
  
  # Determine whether to update the residual variance parameter. Note
  # that the default seting is determined by whether sigma is
  # provided.
  if (missing(update.sigma))
    update.sigma <- update.sigma.default

  # Determine whether to update the mixture variance parameters. By default,
  # these parameters are not updated.
  if (missing(update.sa))
    update.sa <- FALSE
  if (update.sa)
    stop("Estimation of mixture variances not implemented")
  
  # Determine whether to update the mixture weights.
  if (missing(update.w))
    update.w <- TRUE

  # Set initial estimates of variational parameters alpha, ensuring
  # that the smallest value is not less than the "drop threshold" for
  # the mixture components. These parameters are stored as a p x K
  # matrix.
  if (missing(alpha)) {
    alpha <- rand(p,K) + K*drop.threshold
    alpha <- alpha / rep.col(rowSums(alpha),K)
  }
  if (nrow(alpha) != p)
    stop("Input alpha should have as many rows as X has columns")
  if (ncol(alpha) != K)
    stop("Input alpha should have one column for each mixture component")
  if (any(alpha < drop.threshold))
    stop("Initial estimates of \"alpha\" must all be above \"drop.threshold\"")
  
  # Set initial estimates of variational parameters 'mu'. These
  # parameters are stored as a p x K matrix. Note that the first
  # column of this matrix is always zero because it corresponds to the
  # "spike" component.
  if (missing(mu))
    mu <- randn(p,K)
  if (nrow(mu) != p)
    stop("Input mu should have as many rows as X has columns")
  if (ncol(mu) != K)
    stop("Input mu should have one column for each mixture component")
  mu[,1] <- 0

  # (3) PREPROCESSING STEPS
  # -----------------------
  # Adjust the genotypes and phenotypes so that the linear effects of
  # the covariates are removed. This is equivalent to integrating out
  # the regression coefficients corresponding to the covariates with
  # respect to an improper, uniform prior.
  out <- remove.covariate.effects(X,Z,y)
  X   <- out$X
  y   <- out$y
  SZy <- out$SZy
  SZX <- out$SZX
  rm(out)

  # Provide a brief summary of the analysis.
  if (verbose) {
    cat("Fitting variational approximation for linear regression",
        "model with\n")
    cat("mixture-of-normals priors.\n")
    cat(sprintf("samples:      %-6d ",n))
    cat(sprintf("mixture component sd's:         %0.2g..%0.2g\n",
                min(sqrt(sa[2:K])),max(sqrt(sa[2:K]))))
    cat(sprintf("variables:    %-6d ",p))
    cat(sprintf("mixture component drop thresh.: %0.1e\n",drop.threshold))
    cat(sprintf("covariates:   %-6d ",max(0,ncol(Z) - 1)))
    cat(sprintf("fit mixture weights:            %s\n",tf2yn(update.w)))
    cat(sprintf("mixture size: %-6d ",K))
    cat(sprintf("fit residual var. (sigma):      %s\n",tf2yn(update.sigma)))
    cat("intercept:    yes    ")
    cat(sprintf("convergence tolerance      %0.1e\n",tol))
  }

  # Compute a few useful quantities.
  xy <- c(y %*% X)
  d  <- diagsq(X)
  Xr <- c(X %*% rowSums(alpha*mu))
  
  # For each variable and each mixture component, calculate s[i,k],
  # the variance of the regression coefficient conditioned on being
  # drawn from the kth mixture component. Note that first column of
  # "s" is always zero since this corresponds to the "spike" mixture
  # component.
  s <- matrix(0,p,K)
  for (i in 2:K)
    s[,i] <- sigma*sa[i]/(sa[i]*d + 1)

  # Initialize storage for outputs logZ, err and nzw.
  logZ <- rep(0,maxiter)
  err  <- rep(0,maxiter)
  nzw  <- rep(0,maxiter)

  # Initialize the "inactive set"; that is the mixture components with
  # weights that are exactly zero. Also, keep the initial set of
  # mixture variances (sa) and the initial number of mixture
  # components (K). The term "inactive set" is borrowed from "active
  # set methods" in numerical optimization.
  inactive   <- 1:K
  K0         <- K
  w0.penalty <- w.penalty
  sa0        <- sa
  
  # (4) FIT MODEL TO DATA
  # ---------------------
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  if (verbose) {
    cat("       variational    max. --------- hyperparameters ---------\n")
    cat("iter   lower bound  change   sigma  mixture sd's  mix. weights",
        "(drop)\n")
  }
  for (iter in 1:maxiter) {

    # Save the current variational parameters and model parameters.
    alpha0 <- alpha
    mu0    <- mu
    s0     <- s
    sigma0 <- sigma
    w0     <- w

    # (4a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
    logZ0 <- computevarlbmix(Z,Xr,d,y,sigma,sa,w,alpha,mu,s)
    
    # (4b) UPDATE VARIATIONAL APPROXIMATION
    # -------------------------------------
    # Run a forward or backward pass of the coordinate ascent updates.
    if (iter %% 2)
      i <- 1:p
    else
      i <- p:1
    out   <- varbvsmixupdate(X,sigma,sa,w,xy,d,alpha,mu,Xr,i)
    alpha <- out$alpha
    mu    <- out$mu
    Xr    <- out$Xr
    rm(out)

    # (4c) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
    logZ[iter] <- computevarlbmix(Z,Xr,d,y,sigma,sa,w,alpha,mu,s)
    
    # (4d) UPDATE RESIDUAL VARIANCE
    # -----------------------------
    # Compute the approximate maximum likelihood estimate of the residual
    # variable (sigma), if requested. Note that we should also
    # recalculate the variance of the regression coefficients when this
    # parameter is updated. 
    if (update.sigma) {
      sigma <-
        (norm2(y - Xr)^2 + dot(d,betavarmix(alpha,mu,s))
         + sum(colSums(as.matrix(alpha[,-1]*(s[,-1] + mu[,-1]^2)))/sa[-1]))/
           (n + sum(alpha[,-1]))
      for (i in 2:K)
        s[,i] <- sigma*sa[i]/(sa[i]*d + 1)
    }

    # (4e) UPDATE MIXTURE WEIGHTS
    # ---------------------------
    # Compute the approximate penalized maximum likelihood estimate of
    # the mixture weights (w), if requested.
    if (update.w) {
      w <- colSums(alpha) + w.penalty - 1
      w <- w/sum(w)
    }

    # (4f) CHECK CONVERGENCE
    # ----------------------
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum difference
    # between the posterior mixture assignment probabilities at two
    # successive iterations is less than the specified tolerance, or
    # when the variational lower bound has decreased.
    err[iter] <- max(abs(alpha - alpha0))
    nzw[iter] <- K0 - K
    if (verbose) {
      progress.str <-
        sprintf("%04d %+13.6e %0.1e %0.1e %13s [%0.3f,%0.3f] (%d)",
                iter,logZ[iter],err[iter],sigma,
                sprintf("[%0.1g,%0.1g]",sqrt(min(sa[-1])),sqrt(max(sa))),
                min(w),max(w),nzw[iter])
      cat(progress.str)
      cat(rep("\r",nchar(progress.str)))
    }
    if (logZ[iter] < logZ0) {
      logZ[iter] <- logZ0
      err[iter]  <- 0
      sigma      <- sigma0
      w          <- w0
      alpha      <- alpha0
      mu         <- mu0
      s          <- s0
      break
    } else if (err[iter] < tol)
      break

    # (4g) ADJUST INACTIVE SET
    # ------------------------
    # Check if any mixture components should be dropped based on
    # "drop.threshold". Note that the first mixture component (the
    # "spike") should never be droped.
    keep    <- apply(alpha,2,max) >= drop.threshold
    keep[1] <- TRUE
    if (!all(keep)) {

      # At least one of the mixture components satisfy the criterion
      # for being dropped, so adjust the inactive set.
      inactive  <- inactive[keep]
      sa        <- sa[keep]
      w         <- w[keep]
      w0        <- w0[keep]
      w.penalty <- w.penalty[keep]
      alpha     <- alpha[,keep]
      alpha0    <- alpha0[,keep]
      mu        <- mu[,keep]
      mu0       <- mu0[,keep]
      s         <- s[,keep]
      s0        <- s0[,keep]
      K         <- length(inactive)
    }
  }
  if (verbose)
    cat("\n")

  # (6) CREATE FINAL OUTPUT
  # -----------------------
  K   <- K0
  fit <- list(n = n,mu.cov = NULL,update.sigma = update.sigma,
              update.sa = update.sa,update.w = update.w,
              w.penalty = w0.penalty,drop.threshold = drop.threshold,
              sigma = sigma,sa = sa0,w = rep(0,K),alpha = matrix(0,p,K),
              mu = matrix(0,p,K),s = matrix(0,p,K),lfsr = NULL,
              logZ = logZ[1:iter],err = err[1:iter],nzw = nzw[1:iter])
  fit$w[inactive]      <- w
  fit$alpha[,inactive] <- alpha
  fit$mu[,inactive]    <- mu
  fit$s[,inactive]     <- s
  
  # Compute the posterior mean estimate of the regression
  # coefficients for the covariates under the current variational
  # approximation.
  fit$mu.cov <- c(SZy - SZX %*% rowSums(alpha * mu))

  # Compute the local false sign rate (LFSR) for each variable.
  fit$lfsr <- computelfsrmix(alpha,mu,s)
  
  # Add row names to some of the outputs.
  rownames(fit$alpha) <- colnames(X)
  rownames(fit$mu)    <- colnames(X)
  rownames(fit$s)     <- colnames(X)
  names(fit$lfsr)     <- colnames(X)
  
  # Declare the return value as an instance of class 'varbvsmix'.
  class(fit) <- c("varbvsmix","list")
  return(fit)
}

# ----------------------------------------------------------------------
# betavarmix(p,mu,s) returns variances of variables drawn from mixtures of
# normals. Each of the inputs is a n x k matrix, where n is the number of
# variables and k is the number of mixture components. Specifically,
# variable i is drawn from a mixture in which the jth mixture component is
# the univariate normal with mean mu[i,j] and variance s[i,j].
#
# Note that the following two lines should return the same result when k=2
# and the first component is the "spike" density with zero mean and
# variance.
#
#   y1 <- betavar(p,mu,s)
#   y2 <- betavarmix(c(1-p,p),cbind(0,mu),cbind(0,s))
#
betavarmix <- function (p, mu, s)
  rowSums(p*(s + mu^2)) - rowSums(p*mu)^2

# ----------------------------------------------------------------------
# Compute the lower bound to the marginal log-likelihood.
computevarlbmix <- function (Z, Xr, d, y, sigma, sa, w, alpha, mu, s) {

  # Get the number of samples (n), variables (p) and mixture
  # components (K).
  n <- length(y)
  p <- length(d)
  K <- length(w)

  # Compute the variational lower bound.
  out <- (-n/2*log(2*pi*sigma)
          - determinant(crossprod(Z),logarithm = TRUE)$modulus/2
          - (norm2(y - Xr)^2 + dot(d,betavarmix(alpha,mu,s)))/(2*sigma))
  for (i in 1:K)
    out <- (out + sum(alpha[,i]*log(w[i] + eps)) 
                - sum(alpha[,i]*log(alpha[,i] + eps)))
  for (i in 2:K)
    out <- (out + (sum(alpha[,i]) + sum(alpha[,i]*log(s[,i]/(sigma*sa[i]))))/2
                - sum(alpha[,i]*(s[,i] + mu[,i]^2))/(sigma*sa[i])/2)
  return(out)
}

# ----------------------------------------------------------------------
# Compute the local false sign rate (LFSR) for each variable. This
# assumes that the first mixture component is a "spike" (that is, a
# normal density with a variance approaching zero).
computelfsrmix <- function (alpha, mu, s) {

  # Get the number of variables (p) and the number of mixture
  # components (k).
  p <- nrow(alpha)
  k <- ncol(alpha)

  # For each variable, get the posterior probability that the
  # regression coeffiicient is exactly zero.
  p0 <- alpha[,1]

  # For each variable, get the posterior probability that the
  # regression coefficient is negative.
  if (k == 2)
    pn <- alpha[,2] * pnorm(0,mu[,2],sqrt(s[,2]))
  else
    pn <- rowSums(alpha[,-1] * pnorm(0,mu[,-1],sqrt(s[,-1])))
  
  # Compute the local false sign rate (LFSR) following the formula
  # given in the Biostatistics paper, "False discovery rates: a new
  # deal".
  lfsr     <- rep(0,p)
  b        <- pn > 0.5*(1 - p0)
  lfsr[b]  <- 1 - pn[b]
  lfsr[!b] <- p0[!b] + pn[!b]

  return(lfsr)
}
