create.snps <- function (p, n) {
  # Generates minor allele frequencies and additive effects for
  # genetic loci (specifically, these are single nucleotide
  # polymorphisms, or SNPs for short). Additive effects are generated
  # from the standard normal, and minor allele frequencies are uniform
  # between 0.05 and 0.5.
  #
  # Args:
  #   p  Number of SNPs. 
  #   n  Number of causal SNPs.
  #
  # Returns a list containing two components:
  #   maf   Vector with minor allele frequencies of SNPs.
  #   beta  Vector with additive effects of SNPs.
  
  # Generate additive effects for the SNPs, such that N of them have a
  # nonzero effect on the trait.
  S       <- sample(p)
  S       <- S[1:n]
  beta    <- rep(0,p)
  beta[S] <- rnorm(n)

  # Generate the minor allele frequencies. They are uniform on [0.05,0.5].
  maf <- 0.05 + 0.45 * runif(p)

  # Output the minor allele frequencies and additive effects.
  snps <- list(maf = maf,beta = beta)
  return(snps)
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
  data <- list(X = X,y = y)
  return(data)
}

# This is the description of the MATLAB version of this function:
#
# [LNZ,ALPHA,MU,S] = VARBVS(X,Y,SIGMA,SA,LOGODDS) implements the
# fully-factorized variational approximation for Bayesian variable selection
# in linear regression. It finds the "best" fully-factorized variational
# approximation to the posterior distribution of the coefficients in a
# linear regression model of a continuous outcome (quantitiative trait),
# with spike and slab priors on the coefficients. By "best", we mean the
# approximating distribution that locally minimizes the Kullback-Leibler
# divergence between the approximating distribution and the exact posterior.
#
# The required inputs are as follows. Input X is an N x P matrix of
# observations about the variables (or features), where N is the number of
# samples, and P is the number of variables. Y is the vector of observations
# about the outcome; it is a vector of length N. To account for an
# intercept, Y and X must be centered beforehand so that Y and each column
# of X has a mean of zero. 
#
# Note that this routine is implemented with the assumption that the data X
# is single floating-point precision (type HELP SINGLE), as opposed to the
# MATLAB default of double precision. This is useful for large data sets,
# because single precision requires half of the number of bits as double
# floating-point precision. If X is provided in another numerical
# representation, it is immediately converted to SINGLE.
#
# Inputs SIGMA, SA and LOGODDS are the hyperparameters. SIGMA and SA are
# scalars. SIGMA specifies the variance of the residual, and SA*SIGMA is the
# prior variance of the coefficients. LOGODDS is the prior log-odds of
# inclusion for each variable. It is equal to LOGODDS = LOG(Q./(1-Q)), where
# Q is the prior probability that each variable is included in the linear
# model of Y (i.e. the prior probability that its coefficient is not zero).
# LOGODDS may either be a scalar, in which case all the variables have the
# same prior inclusion probability, or it may be a vector of length P.
#
# There are four outputs. Output scalar LNZ is the variational estimate of
# the marginal log-likelihood given the hyperparameters SIGMA, SA and
# LOGODDS.
#
# Outputs ALPHA, MU and S are the parameters of the variational
# approximation and, equivalently, variational estimates of posterior
# quantites: under the variational approximation, the ith regression
# coefficient is normal with probability ALPHA(i); MU(i) and S(i) are the
# mean and variance of the coefficient given that it is included in the
# model. Outputs ALPHA, MU and S are all column vectors of length P.
#
# VARBVS(...,OPTIONS) overrides the default behaviour of the algorithm. Set
# OPTIONS.ALPHA and OPTIONS.MU to override the random initialization of
# variational parameters ALPHA and MU. Set OPTIONS.VERBOSE = FALSE to turn
# off reporting the algorithm's progress.
varbvs <- function (X, y, sigma, sa, logodds, alpha, mu, verbose = TRUE) {

  # Convergence is reached when the maximum relative distance between
  # successive updates of the variational parameters is less than this
  # quantity.
  tolerance <- 1e-4  

  # CHECK INPUTS.
  if (!is.double(X))
    X <- double(X)
  X = single(X)
  
  # Get the number of samples (n) and variables (p).
  n <- nrows(X)
  p <- ncols(X)

  # TAKE CARE OF OPTIONAL INPUTS.
  # Set initial estimates of variational parameters.
  if (is.null(alpha)) {
    alpha <- runif(p)
    alpha <- alpha / sum(alpha)
  }
  if (length(alpha) != p)
    stop("ALPHA and MU must be vectors of length P")
  if (is.null(mu))
    mu <- rnorm(p)    
  if (length(mu) != p)
    stop("MU must be a vector of length P")
}

  # Determine whether to display the algorithm's progress.

##   if isfield(options,'verbose')
##     verbose = options.verbose;
##   else
##     verbose = true;
##   end
##   clear options
  
##   % INITIAL STEPS.
##   % Compute a few useful quantities. Here I calculate X'*Y as (Y'*X)' to
##   % avoid storing the transpose of X, since X may be large.
##   xy = double(y'*X)';
##   d  = diagsq(X);
##   Xr = double(X*(alpha.*mu));
  
##   % Calculate the variance of the coefficients.
##   s = sa*sigma./(sa*d + 1);
  
##   % Repeat until convergence criterion is met.
##   lnZ  = -Inf;
##   iter = 0;
##   if verbose
##     fprintf('       variational    max. incl max.\n');
##     fprintf('iter   lower bound  change vars E[b]\n');
##   end
##   while true

##     % Go to the next iteration.
##     iter = iter + 1;
    
##     % Save the current variational parameters and lower bound.
##     alpha0  = alpha;
##     mu0     = mu;
##     lnZ0    = lnZ;
##     params0 = [ alpha; alpha .* mu ];

##     % UPDATE VARIATIONAL APPROXIMATION.
##     % Run a forward or backward pass of the coordinate ascent updates.
##     if isodd(iter)
##       I = 1:p;
##     else
##       I = p:-1:1;
##     end
##     [alpha mu Xr] = varbvsupdate(X,sigma,sa,logodds,xy,d,alpha,mu,Xr,I);
    
##     % COMPUTE VARIATIONAL LOWER BOUND.
##     % Compute the lower bound to the marginal log-likelihood.
##     lnZ = intlinear(Xr,d,y,sigma,alpha,mu,s) ...
## 	  + intgamma(logodds,alpha) ...
## 	  + intklbeta(alpha,mu,s,sigma*sa);
    
##     % CHECK CONVERGENCE.
##     % Print the status of the algorithm and check the convergence criterion.
##     % Convergence is reached when the maximum relative difference between
##     % the parameters at two successive iterations is less than the specified
##     % tolerance, or when the variational lower bound has decreased. I ignore
##     % parameters that are very small.
##     params = [ alpha; alpha .* mu ];
##     I      = find(abs(params) > 1e-6);
##     err    = relerr(params(I),params0(I));
##     if verbose
##       fprintf('%4d %+13.6e %0.1e %4d %0.2f\n',iter,lnZ,max(err),...
## 	      round(sum(alpha)),max(abs(alpha.*mu)));
##     end
##     if lnZ < lnZ0
##       alpha = alpha0;
##       mu    = mu0;
##       lnZ   = lnZ0;
##       break
##     elseif max(err) < tolerance
##       break
##     end
##   end

# [ALPHA,MU,XR] = VARBVSUPDATE(X,SIGMA,SA,LOGODDS,XY,D,ALPHA0,MU0,XR0,I)
# runs a single iteration of the coordinate ascent updates to maximize the
# variational lower bound for Bayesian variable selection in linear
# regression. It adjusts the fully-factorized variational approximation to
# the posterior distribution of the coefficients in a linear regression
# model of a continuous outcome (quantitiative trait), with spike and slab
# priors on the coefficients.
#
# All inputs to this function are required. Input X is an N x P matrix of
# observations about the variables (or features), where N is the number of
# samples, and P is the number of variables. Input XY = X'*Y, where Y is the
# vector of observations about the outcome. Crucially, to account for an
# intercept, Y and X must be centered beforehand so that Y and each column
# of X has a mean of zero.
#
# This routine is implemented with the assumption that X is a single
# floating-point precision matrix (type HELP SINGLE), as opposed to the
# MATLAB default of double precision. This is useful for large data sets,
# because single precision requires half of the number of bits as double
# floating-point precision. If X is provided in another numerical
# representation, an error is reported.
#
# Inputs SIGMA, SA and LOGODDS specify the hyperparameters. SIGMA and SA are
# scalars. SIGMA specifies the variance of the residual, and SA*SIGMA is the
# prior variance of the regression coefficients. LOGODDS is the prior
# log-odds of inclusion for each variable. It is equal to LOGODDS =
# LOG(Q./(1-Q)), where Q is the prior probability that each variable is
# included in the linear model of Y. LOGODDS is a vector of length P.
#
# Inputs ALPHA0, MU0 are the current parameters of the variational
# approximation; under the variational approximation, the ith regression
# coefficient is normal with probability ALPHA0(i), and MU0(i) is the mean
# of the coefficient given that it is included in the model. Inputs XR0 and
# D must be XR0 = X*(ALPHA0.*MU0) and D = DIAG(X'*X).
#
# Input I specifies the order in which the coordinates are updated. It may
# be a vector of any length. Each entry of I must be an integer between 1
# and P.
#
# There are three outputs. Output vectors ALPHA and MU are the updated
# variational parameters, and XR = X*(ALPHA.*MU). The computational
# complexity of VARBVSUPDATE is O(N*LENGTH(I)).
varbvsupdate <- function (X, sigma, sa, logodds, xy,
                          d, alpha0, mu0, Xr0, I) {
  # *** DESCRIPTION OF FUNCTION GOES HERE ***
  
  # CHECK THE INPUTS.
  # Check input X.
  if (!is.double(X) || !is.matrix(X))
    stop("Input argument X must be a double-precision matrix")
  
  # Get the number of samples (n) and variables (p).
  n <- nrows(X)
  p <- ncols(X)

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

  # Check input I.
  if (sum(I < 1 | I > p) > 0)
      stop("Input I contains invalid variable indices")

  # Execute the C routine, and return the results.
  result <- .C("varbvsupdateR",
               n       = as.integer(n),
               p       = as.integer(p),
               m       = as.integer(length(I)),
               X       = X,
               sigma   = as.double(sigma),
               sa      = as.double(sa),
               logodds = as.double(logodds),
               xy      = as.double(xy),
               d       = as.double(d),
               alpha0  = as.double(alpha0),
               mu0     = as.double(mu0),
               Xr0     = as.double(Xr0),
               I       = as.integer(I),
               alpha   = as.double(alpha0),
               mu      = as.double(mu0),
               Xr      = as.double(Xr0),
               x       = double(n),DUP = FALSE)
  return(list(alpha = result$alpha,
              mu    = result$mu,
              Xr    = result$Xr))
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
