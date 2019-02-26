# This script illustrates how to simulate data from the BSLMM model
# (Zhou, Carbonetto & Stephens, PLoS Genetics, 2013).
source("roots2.R")

# SCRIPT PARAMETERS
# -----------------
n  <- 800   # Number of samples.
p  <- 2000  # Number of variables (these are SNPs).
r  <- 0.5   # Proportion of variance in trait explained by SNPs.
            # This should be a number between 0 and 1.
d  <- 0.3   # Proportion of additive genetic variance due to large effects.
            # This should be a number between 0 and 1.
na <- 20    # Number of "large" (QTL) effects.

# Initialize the sequence of pseudorandom numbers.
set.seed(1)

# GENERATE DATA
# -------------
# Generate the minor allele frequencies so that they are uniform over range
# [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
# (i.e., unlinked), according to the specified minor allele frequencies.
cat("Generating genotype data.\n")
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) +
       (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)

# Center the columns of X.
X <- X - matrix(colMeans(X),n,p,byrow = TRUE)

# Generate the "small" (polygenic) effects for the SNPs.
cat("Generating phenotype data.\n")
u <- randn(p,1)

# Generate the "large" (QTL) additive effects for the SNPs.
i    <- sample(p,na)
b    <- rep(0,p)
b[i] <- rnorm(na)

# Generate random labels for the markers.
colnames(X) <- paste0("rs",sample(1e6,p))

# Adjust the additive effects so that we control for (1) the
# proportion of additive genetic variance that is due to QTL effects
# (d), and (2) the proportion of variance explained (r). That is, we
# adjust b and u so that
#
#   r = a/(a+1)
#   d = c/a,
#
# where I've defined 
#
#   a = (u + b)'*cov(X)*(u + b),
#   c = b'*cov(X)*b.
#
st <- r/(1-r) * d/var(drop(X %*% b))
b  <- sqrt(st) * b
if (d == 0) {

  # All variance is explained by the "small" (polygenic) effects.  
  sa <- r/(1-r)/var(drop(X %*% u))
} else if (d == 1) {

  # All variance is explained by the "large" (QTL) effects.
  sa <- 0
} else {
  sa <- max(roots2(var(drop(X %*% u)),
                   2*sum((X %*% b) * (X %*% u))/n,
                   var(drop(X %*% b)) - r/(1-r)))^2
}
u <- sqrt(sa) * u
    
# Generate the quantitative trait measurements.
y <- drop(X %*% (u + b) + rnorm(n))

# Check that the data were simulated correctly.
cat("Summarizing simulated data:\n")
cat(sprintf(" - Proportion of variance explained by SNPs: %0.3f\n",
            var(drop(X %*% (b + u)))/var(y)))
cat(" - Proportion of additive genetic variance due to large effects: ")
cat(sprintf("%0.3f\n",var(drop(X %*% b))/var(drop(X %*% (b + u)))))
