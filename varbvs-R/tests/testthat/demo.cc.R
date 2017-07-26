# This script illustrates 'varbvs' for genome-wide mapping of a binary
# (e.g., case-control) trait in a simulated data set in which all the
# genetic markers are uncorrelated with each other (i.e., they are
# "unlinked").
#
# Note that variable "covariates" must be specified before running
# this script.
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
n  <- 1400  # Number of samples (subjects).
p  <- 1000  # Number of variables (genetic markers).
na <- 10    # Number of markers that affect the binary outcome.
sa <- 0.2   # Variance of log-odds ratios.
p1 <- 0.5   # Target proportion of subjects that are cases (y = 1).

# Names of covariates. Example:
# 
#   covariates <- c("age","weight")
#

# Candidate values for the prior log-odds of inclusion.
logodds <- seq(-3,-1.5,0.5)

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
X   <- matrix(as.double(X),n,p,byrow = TRUE)

# Generate additive effects for the markers so that exactly na of them
# have a nonzero effect on the trait.
i       <- sample(p,na)
beta    <- rep(0,p)
beta[i] <- sqrt(sa)*rnorm(na)

# Generate random labels for the markers.
colnames(X) <- paste0("rs",sample(1e6,p))

# Generate the covariate data (Z), and the linear effects of the
# covariates (u).
m <- length(covariates)
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
y <- as.double(runif(n) < sigmoid(w))

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a logistic regression model of
# a binary outcome (case-control status), with spike and slab priors
# on the coefficients.
cat("2. FITTING MODEL TO DATA.\n")
fit <- varbvs(X,Z,y,"binomial",logodds = logodds,n0 = 0)
cat("\n")

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("3. SUMMARIZING RESULTS.\n")
print(summary(fit))
cat("\n")

# EVALUATE MODEL PREDICTIONS
# --------------------------
# Compute estimates of the binary trait using the fitted model, and
# compare against the observed values.
cat("4. EVALUATING FITTED MODEL.\n")
y.fit <- predict(fit,X,Z)
cat("Comparison of observed case-control status against estimated outcome:\n")
print(table(y = factor(y),y.fit = factor(y.fit)))
