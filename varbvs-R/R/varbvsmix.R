# Fit linear regression model with mixture-of-normals prior using
# variational approximation. See varbvsmix.Rd for details.
varbvsmix <- function (X, Z, y, sa, sigma, w, alpha, mu, update.sigma,
                       update.sa, update.w, w.penalty, tol = 1e-4,
                       maxiter = 1e4, verbose = TRUE) {
    
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
    w.penalty <- rep(1/K,K)
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

  # Set initial estimates of variational parameters alpha. These
  # parameters are stored as a p x K matrix.
  if (missing(alpha)) {
    alpha <- rand(p,K)
    alpha <- alpha / rep.row(colSums(alpha),p)
  }
  if (nrow(alpha) != p)
    stop("Input alpha should have as many rows as X has columns")
  if (ncol(alpha) != K)
    stop("Input alpha should have one column for each mixture component")

  # Set initial estimates of variational parameters 'mu'. These
  # parameters are stored as a p x K matrix. Note that the first
  # column of this matrix is always zero because it corresponds to the
  # "spike" component.
  if (missing(mu))
    mu <- randn(p,K)
  if (nrow(mu) != p)
    stop("Input mu must have as many rows as X has columns")
  if (ncol(mu) != K)
  mu[,1] <- 0

  # (3) PREPROCESSING STEPS
  # -----------------------
  # Adjust the genotypes and phenotypes so that the linear effects of
  # the covariates are removed. This is equivalent to integrating out
  # the regression coefficients corresponding to the covariates with
  # respect to an improper, uniform prior; see Chipman, George and
  # McCulloch, "The Practical Implementation of Bayesian Model
  # Selection," 2001.
  #
  # Here I compute two quantities that are used here to remove linear
  # effects of the covariates (Z) on X and y, and later on (in
  # function "outerloop"), to efficiently compute estimates of the
  # regression coefficients for the covariates.
  SZy <- solve(crossprod(Z),c(y %*% Z))
  SZX <- solve(crossprod(Z),t(Z) %*% X)
  if (ncol(Z) == 1) {
    X <- X - rep.row(colMeans(X),n)
    y <- y - mean(y)
  } else {

    # The equivalent expressions in MATLAB are  
    #
    #   y = y - Z*((Z'*Z)\(Z'*y))
    #   X = X - Z*((Z'*Z)\(Z'*X))  
    #
    # This should give the same result as centering the columns of X
    # and subtracting the mean from y when we have only one
    # covariate, the intercept.
    y <- y - c(Z %*% SZy)
    X <- X - Z %*% SZX
  }

  # Provide a brief summary of the analysis.
  if (verbose) {
    cat("Fitting variational approximation for linear regression",
        "model with\n")
    cat("mixture-of-normals priors.\n")
    cat(sprintf("samples:      %-6d ",n))
    cat(sprintf("mixture component sd's:    %0.2g..%0.2g\n",
                min(sqrt(sa[2:K])),max(sqrt(sa[2:K]))))
    cat(sprintf("variables:    %-6d ",p))
    cat(sprintf("fit mixture variances:     %s\n",tf2yn(update.sa)))
    cat(sprintf("covariates:   %-6d ",max(0,ncol(Z) - 1)))
    cat(sprintf("fit mixture weights:       %s\n",tf2yn(update.w)))
    cat(sprintf("mixture size: %-6d ",K))
    cat(sprintf("fit residual var. (sigma): %s\n",tf2yn(update.sigma)))
    cat("intercept:    yes    ")
    cat(sprintf("convergence tolerance      %0.1e\n",tol))
  }

  # (4) FIT MODEL TO DATA
  # ---------------------
  # Repeat until convergence criterion is met, or until the maximum
  # number of iterations is reached.
  if (verbose) {
    cat("       variational    max. --------- hyperparameters ---------\n")
    cat("iter   lower bound  change   sigma  mixture sd's  mix. weights\n")
  }
  for (iter in 1:maxiter) {

    # Save the current variational parameters and model parameters.
    alpha0 <- alpha
    mu0    <- mu
    s0     <- s
    sigma0 <- sigma
    sa0    <- sa
    w0     <- w
    
    # (4a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
    logw0 <- computevarlbmix(Z,Xr,d,y,sigma,sa,q,alpha,mu,s)
    
    # (4b) UPDATE VARIATIONAL APPROXIMATION
    # -------------------------------------
    # Run a forward or backward pass of the coordinate ascent updates.
    if (iter %% 2)
      i <- 1:p
    else
      i <- p:1
    # TO DO.

    # (4c) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    # --------------------------------------------
    # Compute the lower bound to the marginal log-likelihood.
    # TO DO.
    
    # (4d) UPDATE RESIDUAL VARIANCE
    # -----------------------------
    # Compute the approximate maximum likelihood estimate of the residual
    # variable (sigma), if requested. Note that we should also
    # recalculate the variance of the regression coefficients when this
    # parameter is updated. 
    if (update_sigma) {
      # TO DO.
    }

    # (4e) UPDATE MIXTURE WEIGHTS
    # ---------------------------
    # Compute the approximate penalized maximum likelihood estimate of
    # the mixture weights (w), if requested.
    if (update.w) {
      # TO DO.
    }

    # (4f) CHECK CONVERGENCE
    # ----------------------
    # Print the status of the algorithm and check the convergence
    # criterion. Convergence is reached when the maximum difference
    # between the posterior inclusion probabilities at two successive
    # iterations is less than the specified tolerance, or when the
    # variational lower bound has decreased.
    # TO DO.
  }
  if (verbose)
    cat("\n")
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
