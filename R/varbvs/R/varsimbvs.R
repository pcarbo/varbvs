varsimbvs <- function (X, y, sigma, sa, log10q, a, b, ca) {
  
  # These two parameters specify the inverse gamma prior on the
  # variance of the residual (sigma). I use the uninformative prior,
  # which is the limit of the inverse gamma as the scale and shape
  # parameters approach zero.
  as <- 0.01
  bs <- 0.01

  # Check input X.
  if (!is.matrix(X))
    stop("Input argument 'X' must be a matrix")

  # Get the number of samples (n), the number of variables (p), and the
  # number of combinations of the hyperparameters (ns).
  n  <- nrow(X)
  p  <- ncol(X)
  ns <- length(sigma)

  # Inputs a, b and c must be scalars.
  if (!is.scalar(a) || !is.scalar(b) || !is.scalar(ca))
    stop("Arguments 'a', 'b' and 'ca' must be scalars")

  # Inputs sigma, sa and log10q must have the same number of elements.
  if (length(sa) != ns || length(log10q) != ns)
    stop("Inputs 'sigma', 'sa' and 'log10q' must be the same size")
  
  # Get the sum of the sample variances.
  sx <- sum(apply(X,2,sd)^2)
  
  # Get the settings for the prior inclusion probabilities.
  q1 <- 10^log10q

  # Initialize storage for the marginal log-likelihoods (lnZ), the
  # log-importance weights (logw), variational estimates of the
  # posterior inclusion probabilities (alpha), and variational
  # estimates of the posterior mean coefficients (mu).
  lnZ     <- sigma
  logw    <- sigma
  lnZ[ ]  <- NA
  logw[ ] <- NA
  alpha   <- matrix(0,p,ns)
  mu      <- matrix(0,p,ns)

  dimnames(alpha) <- list(sample=NULL,importance.sample=NULL)
  dimnames(mu)    <- list(sample=NULL,importance.sample=NULL)
  
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
    result     <- varbvsoptimize(X,y,sigma[i],sa[i],logit(q1[i]),
                                 alpha0,mu0,verbose=FALSE)
    lnZ[i]     <- result$lnZ
    alpha[ ,i] <- result$alpha
    mu[ ,i]    <- result$mu

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
  alpha <- c(alpha %*% c(w))
  mu    <- c(mu    %*% c(w))

  return(list(w=w,alpha=alpha,mu=mu))
}
