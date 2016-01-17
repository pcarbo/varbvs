# This script illustrates 'varbvs' for genome-wide mapping of a binary
# (e.g., case-control) trait in a simulated data set in which all the
# genetic markers are uncorrelated with each other (i.e., they are
# "unlinked").
#
# NOTE: Move this file later to the "demo" directory.
#
source("misc.R")
source("varbvs.R")
source("varbvsbin.R")
source("varbvsbinupdate.R")
dyn.load("../src/diagsqr.so")
dyn.load("../src/varbvsr.so")

# SCRIPT PARAMETERS
# -----------------
n  <- 1500  # Number of samples (subjects).
p  <- 2000  # Number of variables (genetic markers).
m  <- 2     # Number of covariates (m >= 0).
na <- 20    # Number of markers that affect the binary outcome.
sa <- 0.15  # Variance of log-odds ratios.
p1 <- 0.25  # Target proportion of subjects that are cases (y = 1).

# Names of covariates.
covariates <- c("age","weight")

# Candidate values for the prior log-odds of inclusion.
logodds <- seq(-3,-1,0.1)

# Set the random number generator seed.
set.seed(1)

# GENERATE DATA SET
# -----------------
# Generate the minor allele frequencies so that they are uniform over
# range [0.05,0.5]. Then simulate genotypes assuming all markers are
# uncorrelated (i.e., unlinked), according to the specified minor
# allele frequencies.
cat("1. GENERATING DATA SET.\n")
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) +
       (runif(n*p) < maf)
X   <- matrix(X,n,p,byrow = TRUE)
storage.mode(X) <- "double"

# Generate additive effects for the markers so that exactly na of them
# have a nonzero effect on the trait.
i       <- sample(p,na)
beta    <- rep(0,p)
beta[i] <- rnorm(na)

# Generate random labels for the markers.
colnames(X) <- paste0("rs",sample(1e6,p))

# Generate the covariate data (Z), and the linear effects of the
# covariates (u).
if (m > 0) {
  Z <- randn(n,m)
  u <- rnorm(m)
  colnames(Z) <- covariates
} else {
  Z <- NULL
}

# For each sample, calculate the probability of being a case (y = 1).
mu <- logit(p1)
w  <- mu + X %*% beta
if (m > 0)
  w <- w + Z %*% u

# Simulate the binary trait (case-control status) as a coin toss with
# success rates given by the logistic regression.
y <- runif(n) < sigmoid(w)
y <- as.double(y)

# TEMPORARY.
out <- varbvsbin(X,y,1,rep(log(na/p),p),runif(p),rnorm(p),rep(1,n))

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a logistic regression model of
# a binary outcome (case-control status), with spike and slab priors
# on the coefficients.
# TO DO.
