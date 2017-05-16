# This script illustrates the "varbvsmix" function on a simulated data
# set, in which all candidate variables (predictors) are uncorrelated.
library(varbvs)  

# SCRIPT PARAMETERS
# -----------------
n  <- 1000  # Number of samples.
p  <- 2000  # Number of variables (genetic markers).
se <- 4     # Variance of residual.

# Names of covariates.
covariates <- c("age","weight")

# The standard deviations and mixture weights used to simulate the
# additive effects on the quantitative trait. Note that the first
# mixture component must have a standard deviation of exactly zero.
sd <- c(0,   0.1, 0.2, 0.5)
w  <- c(0.95,0.03,0.01,0.01)

# Set the random number generator seed.
set.seed(1)

# GENERATE DATA SET
# -----------------
# Get the number of covariates.
m <- length(covariates)

# Generate the minor allele frequencies so that they are uniform over
# range [0.05,0.5]. Then simulate genotypes assuming all markers are
# uncorrelated (i.e., unlinked), according to the specified minor
# allele frequencies.
cat("1. GENERATING DATA SET.\n")
cat("Data simulation settings:\n")
cat(sprintf("  - Num. data samples       %d\n",n))
cat(sprintf("  - Num. covariates         %d\n",m))
cat(sprintf("  - Num. variables (SNPs)   %d\n",p))
cat(sprintf("  - Num. mixture components %d\n",length(w)))
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) +
       (runif(n*p) < maf)
X   <- matrix(as.double(X),n,p,byrow = TRUE)

# Generate additive effects according to the specified standard
# deviations (sd) and mixture weights (w).
k    <- sample(length(w),p,replace = TRUE,prob = w)
beta <- sd[k] * rnorm(p)

# Generate random labels for the markers.
colnames(X) <- paste0("rs",sample(1e6,p))

# Generate a random intercept.
mu <- rnorm(1)

# Generate the covariate data (Z), and the linear effects of the
# covariates (u).
if (m > 0) {
  Z <- randn(n,m)
  u <- rnorm(m)
  colnames(Z) <- covariates
} else {
  Z <- NULL
}

# Generate the quantitative trait measurements.
y <- mu + X %*% beta + sqrt(se)*rnorm(n)
if (m > 0)
  y <- y + Z %*% u
y <- c(y)

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of
# the quantitative trait (Y), with the mixture-of-normals prior on the
# coefficients.
cat("2. FITTING MODEL TO DATA.\n")
fit <- varbvsmix(X,Z,y,sd^2)

