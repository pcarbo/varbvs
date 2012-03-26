create.snps <- function (p, n) {
  # Generates minor allele frequencies and additive effects for
  # genetic loci (specifically, these are single nucleotide
  # polymorphisms, or SNPs for short). Additive effects are generated
  # from the standard normal, and minor allele frequencies are uniform
  # between 0.05 and 0.5.

  # Check the inputs.
  if (!is.scalar(p))
    stop("Invalid 'p' argument")
  if (!is.scalar(n))
    stop("Invalid 'n' argument")
    
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

create.data <- function (maf, beta, sigma, n) {
  # Generates samples of the genotypes and quantitative trait (the
  # continuous outcome Y) according to the specified SNP minor allele
  # frequencies and additive effects. Genotypes are generated from a
  # binomial distribution with success rates given by the minor allele
  # frequencies. Observations about the quantitative trait are
  # generated according to y = X*beta + e, where the residual e is
  # normal with mean zero and covariance sigma*I.

  # Check the inputs.
  maf  <- c(maf)
  beta <- c(beta)
  if (!is.numeric(maf))
    stop("Invalid 'maf' argument")
  if (!is.numeric(beta))
    stop("Invalid 'beta' argument")
  if (!is.numeric(sigma))
    stop("Invalid 'sigma' argument")
  if (!is.numeric(n))
    stop("Invalid 'n' argument")
  if (length(maf) != length(beta))
    stop("Vectors 'maf' and 'beta' must be the same length")

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
  dimnames(X) <- list(sample=NULL,variable=NULL)
  return(list(X = X,y = y))
}  

