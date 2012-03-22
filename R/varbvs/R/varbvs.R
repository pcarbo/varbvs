# Shorthand constant for machine precision.
eps <- .Machine$double.eps

create.snps <- function (p, n) {
  # Generates minor allele frequencies and additive effects for
  # genetic loci (specifically, these are single nucleotide
  # polymorphisms, or SNPs for short). Additive effects are generated
  # from the standard normal, and minor allele frequencies are uniform
  # between 0.05 and 0.5.
  
  # Generate additive effects for the SNPs, such that N of them have a
  # nonzero effect on the trait.
  S       <- sample(p)
  S       <- S[1:n]
  beta    <- rep(0,p)
  beta[S] <- rnorm(n)

  # Generate the minor allele frequencies. They are uniform on [0.05,0.5].
  maf <- 0.05 + 0.45 * runif(p)

  # Output the minor allele frequencies and additive effects.
  return(list(maf = maf,beta = beta))
}

create.data <- function (snps, sigma, n) {
  # Generates samples of the genotypes and quantitative trait (the
  # continuous outcome Y) according to the specified SNP minor allele
  # frequencies and additive effects. Genotypes are generated from a
  # binomial distribution with success rates given by the minor allele
  # frequencies. Observations about the quantitative trait are
  # generated according to y = X*beta + e, where the residual e is
  # normal with mean zero and covariance sigma*I.
  #
  # Args:
  #   snps$maf   Minor allele frequencies of the SNPs.
  #   snps$beta  Additive effects of the SNPs.
  #   sigma      Variance of residual.
  #   n          Number of samples to generate.
  #
  # Returns a list containing two components:
  #   X  Matrix of genotype data with centered columns.
  #   y  Vector of quantitative trait data (centered so that mean is 0).

  # The the information about the SNPs.
  maf  <- snps$maf
  beta <- snps$beta
  
  # Get the number of SNPs.
  p <- length(maf)
  
  # Simulate genotype data X from an idealized population, according to the
  # specified minor allele frequencies.
  X <- (rnorm(n*p) < rep(maf,n)) +
       (rnorm(n*p) < rep(maf,n))
  X <- matrix(X,n,p)

  # Center the columns of X.
  X <- center.columns(X)

  # Generate the quantitative trait measurements.
  y <- X %*% beta + sqrt(sigma) * rnorm(n)
  y <- c(y)
  
  # Take into account an intercept by centering the outcomes Y to have
  # mean zero.
  y <- y - mean(y)

  # Output the genotype and quantitative trait samples.
  return(list(X = X,y = y))
}

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

varbvsoptimize <- function (X, y, sigma, sa, logodds, alpha0 = NULL,
                            mu0 = NULL, verbose = TRUE) {
  # Implements the fully-factorized variational approximation for
  # Bayesian variable selection in linear regression. It finds the
  # "best" fully-factorized variational approximation to the posterior
  # distribution of the coefficients in a linear regression model of a
  # continuous outcome (quantitiative trait), with spike and slab
  # priors on the coefficients. By "best", we mean the approximating
  # distribution that locally minimizes the Kullback-Leibler
  # divergence between the approximating distribution and the exact
  # posterior.
  #
  # Args:
  #   X        n x p matrix of observations about the variables (or
  #            features), where n is the number of samples
  #            (observations), and p is the number of variables.
  #   y        Vector of observations about the outcome. It is a
  #            vector of length n.
  #   sigma    Variance of the residual.
  #   sa       sa*sigma is the prior variance of the coefficients.
  #   logodds  Prior log-odds of inclusion for each variable. It is
  #            equal to logodds = log(q/(1-q)), where q is the prior
  #            probability that each variable is included in the
  #            linear model of Y (i.e. the prior probability that
  #            its coefficient is not zero). It may either be a
  #            scalar, in which case all the variables have the
  #            same prior inclusion probability, or it may be a
  #            vector of length p.
  #   alpha0   Initial variational estimate of posterior inclusion
  #            probabilities. If alpha0 = NULL, the variational
  #            parameters are initialized at random.
  #   mu0      Initial variational estimate of posterior mean
  #            coefficients. If mu0 = NULL, the variational
  #            parameters are randomly initialized.
  #   verbose  Set verbose = FALSE to turn off reporting the
  #            algorithm's progress. 
  #
  # Returns a list containing four components: ‘alpha’, ‘mu’ and ‘s’
  # are the parameters of the variational approximation and,
  # equivalently, variational estimates of posterior quantites: under
  # the variational approximation, the ith regression coefficient is
  # normal with probability alpha[i]; mu[i] and s[i] are the mean and
  # variance of the coefficient given that it is included in the
  # model. Outputs alpha, mu and s are all column vectors of length p.
  #
  # Output ‘lnZ’ is a scalar giving the variational estimate of
  # the marginal log-likelihood.  
  #
  # Note: to account for an intercept, y and X must be centered
  # beforehand so that vector y and each column of X has a mean of
  # zero.

  # Convergence is reached when the maximum relative distance between
  # successive updates of the variational parameters is less than this
  # quantity.
  tolerance <- 1e-4  

  # CHECK INPUTS.
  # Check input X.
  if (!is.matrix(X))
    stop("Input argument X must be a matrix")
  if (!is.double(X))
    X <- double(X)

  # Get the number of samples (n) and variables (p).
  n <- nrow(X)
  p <- ncol(X)

  # Check input y.
  if (length(y) != n)
    stop("Data X and y do not match")

  # Check inputs sigma and sa.
  if (!is.scalar(sigma) || !is.scalar(sa))
    stop("Input arguments sigma and sa must be scalars")

  # Check input logodds.
  if (is.scalar(logodds))
    logodds <- rep(logodds,p)
  if (length(logodds) != p)
    stop("Input logodds must be a scalar or a vector of length p")

  # TAKE CARE OF OPTIONAL INPUTS.
  # Set initial estimates of variational parameters.
  if (is.null(alpha0)) {
    alpha <- runif(p)
    alpha <- alpha / sum(alpha)
  }
  else
    alpha <- alpha0
  if (is.null(mu0))
    mu <- rnorm(p)
  else
    mu <- mu0
  if (length(alpha) != p || length(mu) != p)
    stop("alpha0 and mu0 must be vectors of length p")
  
  # INITIAL STEPS.
  # Compute a few useful quantities.
  xy <- c(y %*% X)           # xy = X'*y.
  d  <- diagsq(X)            # d  = diag(X'*Y).
  Xr <- c(X %*% (alpha*mu))  # Xr = X*(alpha*mu).

  # Calculate the variance of the coefficients.
  s <- sa*sigma/(sa*d + 1)

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
    params0 <- c(alpha,alpha*mu)

    # UPDATE VARIATIONAL APPROXIMATION.
    # Run a forward or backward pass of the coordinate ascent updates.
    if (is.odd(iter))
      S <- seq(1,p)
    else
      S <- seq(p,1,-1)
    result <- varbvsupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,S)
    alpha  <- result$alpha
    mu     <- result$mu
    Xr     <- result$Xr

    # COMPUTE VARIATIONAL LOWER BOUND.
    # Compute the lower bound to the marginal log-likelihood.
    lnZ <- intlinear(Xr,d,y,sigma,alpha,mu,s) +
           intgamma(logodds,alpha) +
           intklbeta(alpha,mu,s,sigma*sa)

    # CHECK CONVERGENCE.
    # Print the status of the algorithm and check the convergence
    # criterion.  Convergence is reached when the maximum relative
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
      lnZ   <- lnZ0
      break
    }
    else if (max(err) < tolerance)
      break
  }

  # Return the variational estimates.
  return(list(alpha=alpha,mu=mu,s=s,lnZ=lnZ))
}
  
varbvsupdate <- function (X, sigma, sa, logodds, xy, d, alpha0, mu0, Xr0, S) {
  # Runs a single iteration of the coordinate ascent updates to
  # maximize the variational lower bound for Bayesian variable
  # selection in linear regression. It adjusts the fully-factorized
  # variational approximation to the posterior distribution of the
  # coefficients in a linear regression model of a continuous outcome
  # (quantitiative trait), with spike and slab priors on the
  # coefficients.
  #
  # Args:
  #   X        n x p matrix of observations about the variables (or
  #            features), where n is the number of samples, and p is
  #            the number of variables. X must be a double-precision
  #            matrix. (Unlike MATLAB, there are no single-precision
  #            floating-point numbers in R.) 
  #   sigma    Variance of the residual (a scalar).
  #   sa       sa*sigma is the prior variance of the regression 
  #            coefficients (also a scalar).
  #   logodds  Prior log-odds of inclusion for each variable. It is equal
  #            to logodds = log(q./(1-q)), where q is the prior
  #            probability that each variable is included in the linear
  #            model of Y. It is a vector of length p.
  #   xy       Equal to X'*y, where y is the vector of observations
  #            about the outcome.
  #   d        Equal to diag(X'*X).
  #   alpha0   Current variational estimates of the posterior inclusion
  #            probabilities. It is a vector of length p.
  #   mu0      Current variational estimates of the posterior mean
  #            coefficients. It is a vector of length p.
  #   Xr0      Equal to X*(alpha0*mu0).
  #   S        Order in which the coordinates are updated. It is a
  #            vector of any length. Each entry of S must be an
  #            integer between 1 and p.
  #
  # Returns a list containing three components:
  #   alpha    Updated variational estimates of the posterior inclusion
  #            probabilities. 
  #   mu       Updated variational estimates of the posterior mean
  #            coefficients.
  #   Xr       Equal to X*(alpha*mu).
  #
  # Note: to account for an intercept, y and X must be centered
  # beforehand so that y and each column of X has a mean of zero. The
  # computational complexity of this algorithm is O(n*length(S)).
  
  # CHECK THE INPUTS.
  # Check input X.
  if (!is.double(X) || !is.matrix(X))
    stop("Input argument X must be a double-precision matrix")
  
  # Get the number of samples (n), the number of variables (p), and
  # the number of updates to execute (m).
  n <- nrow(X)
  p <- ncol(X)
  m <- length(S)

  # Check inputs sigma and sa.
  if (!is.scalar(sigma) || !is.scalar(sa))
    stop("Input arguments sigma and sa must be scalars")

  # Check input logodds.
  if (is.scalar(logodds))
    logodds <- rep(logodds,p)
  if (length(logodds) != p)
    stop("Input logodds must be a scalar or a vector of length p")

  # Check inputs xy and d.
  if (length(xy) != p || length(d) != p)
    stop("Inputs xy and d must be vectors of length p")

  # Check inputs alpha0 and mu0.
  if (length(alpha0) != p || length(mu0) != p)
    stop("Inputs alpha0 and mu0 must be vectors of length p")

  # Check input Xr0.
  if (length(Xr0) != n)
    stop("Input Xr0 must be a vector of length n")

  # Check input S.
  if (sum(S < 1 | S > p) > 0)
      stop("Input S contains invalid variable indices")

  # Execute the C routine, and return the results in a list object.
  # The only components of the list that change are alpha, mu and Xr.
  # We need to subtract 1 from the indices because R vector start at
  # 1, but C arrays start at 0. Note that I do not attempt to coerce X
  # here; if X is large, it could use a lot of memory to duplicate
  # this matrix. For this same reason, I set DUP = FALSE so that the
  # input arguments are not duplicated.
  result <- .C("varbvsupdateR",
               n       = as.integer(n),      # Number of samples.
               m       = as.integer(m),      # Number of updates.
               X       = X,                  # Matrix of samples.
               sigma   = as.double(sigma),   # Residual variance.
               sa      = as.double(sa),      # Prior variance of coefficients.
               logodds = as.double(logodds), # Prior log-odds.
               xy      = as.double(xy),      # xy = X'*y.
               d       = as.double(d),       # d = diag(X'*X).
               alpha   = as.double(alpha0),  # Posterior inclusion prob's.
               mu      = as.double(mu0),     # Posterior mean coefficients.
               Xr      = as.double(Xr0),     # Xr = X*r.
               S       = as.integer(S-1),    # Updates to perform.
               DUP     = FALSE)
  return(list(alpha = result$alpha,
              mu    = result$mu,
              Xr    = result$Xr))
}

# NORMALIZELOGWEIGHTS(LOGW) takes as input an array of unnormalized
# log-importance weights LOGW and returns normalized importance weights such
# that the sum of the normalized importance weights is equal to one.
normalizelogweights <- function (logw) {

  # We guard against underflow or overflow by adjusting the log-importance
  # weights so that the largest importance weight is one.
  c <- max(c(logw))
  w <- exp(logw - c)

  # Normalize the importance weights.
  w <- w / sum(c(w))
}
  
  # LOGINVGAMMA(X,A,B) returns the logarithm of the probability density
# function of the inverse gamma distribution at elements of X, with
# shape A and scale B.
loginvgamma <- function (x, a, b) {
  f <- a*log(b) - lgamma(a) - (a+1)*log(x+eps) - b/(x+eps)
  return(f)
}
  
# LOGPVE(A,X) returns the logarithm of the density function for X,
# given that R = A*X/(A*X + 1) is uniform on [0,1]. Input A must be a
# positive scalar. This is useful for calculating the prior
# probability of the variance parameter X, when we have a uniform
# prior on R, the proportion of variance explained.
logpve <- function (a, x) {
  r <- a*x/(a*x+1)          # Proportion of variance explained.
  f <- 2*log(r/x) - log(a)  # Log-density.
  return(f)
}
  
intlinear <- function (Xr, d, y, sigma, alpha, mu, s) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the expectation
  # of the linear regression log-likelihood taken with respect to the
  # variational approximation. Inputs Xr and d must be equal to Xr =
  # X*(alpha*mu) and d = diag(X'*X). For a description of the
  # remaining inputs, see function ‘varbvsoptimize’.
  n <- length(y)
  f <- -n/2*log(2*pi*sigma) - norm2(y - Xr)^2/(2*sigma) -
        dot(d,betavar(alpha,mu,s))/(2*sigma)
  return(f)
}
                      
intklbeta <- function (alpha, mu, s, sa) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the negative
  # Kullback-Leibler divergence between the approximating distribution
  # and the prior of the coefficients. Input ‘sa’ specifies the prior
  # variance of the coefficients. (It is not the same as the ‘sa’ used
  # as input to ‘varbvsoptimize’.) See function ‘varbvsoptimize’ for
  # details on the inputs to this function.
  f <- (sum(alpha) + dot(alpha,log(s/sa)) - dot(alpha,s + mu^2)/sa)/2 -
        dot(alpha,log(alpha + eps)) -
        dot(1 - alpha,log(1 - alpha + eps))
  return(f)
}

intgamma <- function (logodds, alpha) {
  # Computes an integral that appears in the variational lower bound
  # of the marginal log-likelihood. This integral is the expectation
  # on the prior inclusion probabilities taken with respect to the
  # variational approximation.
  #
  # Args:
  #   logodds  Scalar, or a vector, specifying the prior log-odds.
  #   alpha    Mixture weights for the variational approximation.
  #
  # See function ‘varbvsoptimize’ for details on the input arguments.

  # This is the same as 
  #
  #    sum(alpha*log(q) + (1-alpha)*log(1-q)).
  #  
  return(sum((alpha-1) * logodds + logsigmoid(logodds)))
}

betavar <- function (p, mu, s) {
  # betavar(p,mu,s) returns the variance of random vector X, in which
  # X[i] is drawn from the normal distribution with probability p[i],
  # and X[i] is zero with probability 1-p[i]. Inputs mu and s specify
  # the mean and variance of the normal density. Inputs p, mu and s
  # must be vectors of the same length. This function is useful for
  # calculating the variance of the coefficients under the
  # fully-factorized variational approximation.

  # Note that this is the same as 
  # 
  #    v = p*(s + mu^2) - (p*mu)^2.
  #
  return(p*(s + (1 - p)*mu^2))
}

dot <- function (x,y) {
  # Return the dot product of vectors x and y.
  return(sum(x*y))
}

norm2 <- function (x) {
  # Return the quadratic norm (2-norm) of vector x.
  return(sqrt(dot(x,x)))
}

logsigmoid <- function (x) {
  # Returns the logarithm of the sigmoid. Use this instead of
  # log(sigmoid(x)) to avoid loss of numerical precision.
  return(-logpexp(-x))
}

logpexp <- function (x) {
  # logpexp(x) returns log(1 + exp(x)). For large x, logpexp(x) should
  # be approximately x. The computation is performed in a numerically
  # stable manner.

  # For large entries, log(1 + exp(x)) is effectively the same as x.
  y <- x

  # Find entries of x that are not large. For these entries, compute
  # log(1 + exp(x)).
  S    <- which(x < 16)
  y[S] <- log(1 + exp(x[S]))
  return(y)
}

logit <- function (x) {
  # The logit function, which is the reverse of the sigmoid function.
  return(log((x + eps)/((1 - x) + eps)))
}
  
diagsq <- function (X, a = NULL) {
  # diagsq(X) returns diag(X'*X).
  # diagsq(X,a) returns diag(X'*diag(a)*X).
  if (is.null(a)) {
    n <- nrow(X)
    a <- rep(1,n)
  }
  y <- c(a %*% X^2)  # y = X^2*a.
  return(y)
}

relerr <- function (x1, x2) {
  # Returns the absolute relative error.
  if (length(x1) == 0 || length(x2) == 0)
    y <- 0
  else
    y <- abs(x1 - x2) / (abs(x1) + abs(x2) + eps)
  return(y)
}

repmat <- function (A,m,n) {
  # Does the same thing as REPMAT(A,m,n) in MATLAB.
  return(kronecker(matrix(1,m,n),A))
}

center.columns <- function (X) {
  # Centers the columns of X so that the entries in each column of X
  # add up to zero.
  mu <- matrix(colMeans(X),1,ncol(X))
  mu <- repmat(mu,nrow(X),1)
  X  <- X - mu
  return(X)
}

grid3d <- function (x,y,z) {
  # Does the same thing as NDGRID(x,y,z) in MATLAB.

  # Get the number of entries in each of the inputs.
  nx <- length(x)
  ny <- length(y)
  nz <- length(z)

  # Initialize the 3-d outputs.
  X <- array(dim=c(nx,ny,nz))
  Y <- X
  Z <- X

  # Set the entries of the 3-d arrays.
  for (i in 1:nx)
    for (j in 1:ny)
      for (k in 1:nz) {
        X[i,j,k] <- x[i]
        Y[i,j,k] <- y[j]
        Z[i,j,k] <- z[k]
      }
  
  grid <- list(X = X,Y = Y,Z = Z)
  return(grid)
}

is.scalar <- function (x) {
  # Returns TRUE if and only if x is a scalar.
  return(length(x) == 1)
}

is.odd <- function (x) {
  # Returns TRUE if and only if x is odd.
  return(x %% 2)
}

