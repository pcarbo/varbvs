## % [W,ALPHA,MU] = VARSIMBVS(X,Y,SIGMA,SA,LOG10Q,A,B,C) runs the full
## % inference procedure for Bayesian variable selection in linear
## % regression. This is a special implementation of the variational inference
## % procedure used in the two simulation studies for the Bayesian Analysis
## % paper. The main distinguishing feature of this procedure is the choice of
## % priors for the hyperparameters of the variable selection model. In
## % addition, we also avoid erratic behaviour in the approximation by first
## % searching for a good initialization of the variational parameters. This
## % inference procedure involves an inner loop and an outer loop. The inner
## % loop consists of running a coordinate ascent algorithm to tighten the
## % variational lower bound given a setting of the hyperparameters (this inner
## % loop is implemented in function VARBVSOPTIMIZE). The outer loop computes
## % importance weights for all combinations of the hyperparameters.
## %
## % Input X is an N x P matrix of observations about the variables (or
## % features), where N is the number of samples, and P is the number of
## % variables. Y is the vector of observations about the outcome; it is a
## % vector of length N. To account for an intercept, Y and X must be centered
## % beforehand so that Y and each column of X has a mean of zero.
## %
## % Note that this routine is implemented with the assumption that the data X
## % is single floating-point precision (type HELP SINGLE), as opposed to the
## % MATLAB default of double precision. This is useful for large data sets,
## % because single precision requires half of the number of bits as double
## % floating-point precision. If X is provided in another numerical
## % representation, it is immediately converted to SINGLE.
## %
## % Inputs SIGMA, SA and LOG10Q specify the hyperparameter settings. These
## % inputs must be arrays of the same size. For each combination of the
## % hyperparameters, we compute an importance weight, and store the result in
## % output W. SIGMA is the residual variance, SIGMA.*SA is the prior variance
## % of the regression coefficients, and LOG10Q is the (base 10) logarithm of
## % the prior inclusion probability.
## %
## % Inputs A, B and C are positive scalars. A and B are the prior sample sizes
## % for the Beta prior on the prior inclusion probability. We assume a uniform
## % prior on the proportion of variance explained, except that we replace the
## % prior inclusion probability the proportion of variance explained by a
## % constant, C. This is done purely for convenience, so that hyperparameter
## % SA does not depend on the prior inclusion probability a priori, making it
## % easier to implement the Markov chain Monte Carlo (MCMC) method. We assume
## % the standard noninformative prior on the residual variance SIGMA.
## %
## % Outputs ALPHA and MU are variational estimates of the posterior inclusion
## % probabilities and posterior mean of the coefficients (given that the
## % variable is included in the model). These variational estimates are
## % averaged over the settings of the hyperparameters.
varsimbvs <- function (X, y, sigma, sa, log10q, a, b, ca) {
  
  # These two parameters specify the inverse gamma prior on the
  # variance of the residual (sigma). I use the uninformative prior,
  # which is the limit of the inverse gamma as the scale and shape
  # parameters approach zero.
  as <- 0.01
  bs <- 0.01

  # Check input X.
  if (!is.matrix(X))
    stop("Input argument X must be a matrix")

  # Get the number of samples (n), the number of variables (p), and the
  # number of combinations of the hyperparameters (ns).
  n  <- nrow(X)
  p  <- ncol(X)
  ns <- length(sigma)

  # Inputs a, b and c must be scalars.
  if (!is.scalar(a) || !is.scalar(b) || !is.scalar(c))
    stop("Input arguments a, b and c must be scalars")

  # Inputs sigma, sa and log10q must have the same number of elements.
  if (length(sa) != ns || length(log10q) != ns)
    stop("Inputs sigma, sa and log10q must be the same size")
  
  # Get the sum of the sample variances.
  sx <- sum(sd(X)^2)
  
  # Get the settings for the prior inclusion probabilities.
  q1 <- 10^log10q

  # Initialize storage for the marginal log-likelihoods (lnZ), the
  # log-importance weights (logw), variational estimates of the
  # posterior inclusion probabilities (alpha), and variational
  # estimates of the posterior mean coefficients (mu).
  lnZ   <- rep(0,ns)
  logw  <- rep(0,ns)
  alpha <- matrix(0,p,ns)
  mu    <- matrix(0,p,ns)
  
  # First get the best initialization for the variational
  # parameters. Repeat for each combination of the hyperparameters.
  cat(sprintf("Finding best initialization for %d combinations ",ns))
  cat("of hyperparameters.\n")
  for (i in 1:ns) {
    cat(sprintf("(%03d) sigma = %4.1f, sa = %0.3f, q = %0.2e",
                i,sigma[i],sa[i],q1[i]))
    cat(rep("\b",1,44),sep="");
    
    # Randomly initialize the variational parameters.
    alpha0 <- runif(p)
    alpha0 <- alpha0 / sum(alpha0)
    mu0    <- rnorm(p)

    # Run the coordinate ascent algorithm.
    result     <- varbvsoptimize(X,y,sigma[i],sa[i],logit(q1[i]),
                                 alpha0,mu0,verbose=FALSE)
    lnZ[i]     <- result$lnZ
    alpha[ ,i] <- result$alpha
    mu[ ,i]    <- result$mu
  }
  cat("\n")
  
  # Choose an initialization common to all the runs of the coordinate ascent
  # algorithm. This is chosen from the hyperparameters with the highest
  # marginal likelihood.
  i      <- which.max(lnZ)
  alpha0 <- c(alpha[ ,i])
  mu0    <- c(mu[ ,i])
  
  # Repeat for each combination of the hyperparameters.
  cat(sprintf("Computing importance weights for %d combinations ",ns))
  cat("of hyperparameters.\n")
  for (i in 1:ns) {
    cat(sprintf("(%03d) sigma = %4.1f, sa = %0.3f, q = %0.2e",
                i,sigma[i],sa[i],q1[i]))
    cat(rep("\b",1,44),sep="");
  
    # Run the coordinate ascent algorithm.
    result <- varbvsoptimize(X,y,sigma[i],sa[i],logit(q1[i]),alpha0,mu0,
                             verbose=FALSE)

    # Compute the log-importance weight. Note that if X ~ p(x), then
    # the probability density of Y = log(X) is proportional to
    # x*p(x). This is useful for calculating the prior for the
    # logarithm of the prior inclusion probability, since we want to
    # calculate the posterior distribution for log10(q), not q.
    logw[i] <- lnZ[i] +                           # Marginal likelihood.
               loginvgamma(sigma[i],as/2,bs/2) +  # Prior on sigma.
               logpve(ca*sx,sa[i]) +              # Prior on sa.
               dbeta(q1[i],a+1,b,log=TRUE)        # Prior on log10q.
  }
  cat("\n")
  
  # Compute the normalized importance weights.
  w <- normalizelogweights(logw)

  # Compute the weighted averages.
  alpha <- c(alpha %*% w)
  mu    <- c(mu    %*% w)

  return(list(w=w,alpha=alpha,mu=mu))
}
