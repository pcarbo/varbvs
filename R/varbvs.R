# This is the description of the MATLAB version of this function:
#
# [MAF,BETA] = CREATESNPS(P,N) generates minor allele frequencies (MAF) and
# additive effects (BETA) for P genetic loci (specifically, these are single
# nucleotide polymorphisms, or SNPs for short). Additive effects are
# generated from the standard normal. Minor allele frequencies are uniform
# between 0.05 and 0.5. Input argument N specifies the number of causal
# SNPs (SNPs that have a nonzero additive effect on the trait).
create.snps <- function (p, n) {
  
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

# This is the description of the MATLAB version of this function:
#
# [X,Y] = CREATEDATA(MAF,BETA,SIGMA,N) generates N samples of the genotypes
# and the quantitative trait (the continuous outcome Y), according to SNP
# minor allele frequencies MAF, and additive effects BETA. Inputs MAF and
# BETA are vectors of length P. The genotype data X is an N x P matrix, and
# the quantitative trait data Y is a column vector of length N. Both X and Y
# are centered so that Y and each column of X has a mean of zero.
#
# Genotypes are generated from a binomial distribution with success rates
# given by the minor allele frequencies. Observations about the quantitative
# trait are generated according to Y = X*BETA + E, where the residual E is
# normal with mean zero and covariance sigma*I.
create.data <- function (snps, sigma, n) {

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

# This replicates the REPMAT function from MATLAB.
repmat <- function (X,m,n) { return(kronecker(matrix(1,m,n),X)) }

# Center the columns of matrix X so that the entries in each column of
# X add up to zero.
center.columns <- function (X) {
  mu <- matrix(colMeans(X),1,ncol(X))
  mu <- repmat(mu,nrow(X),1)
  X  <- X - mu
  return(X)
}
