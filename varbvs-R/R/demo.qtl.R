# This script illustrates usage of the function varbvs for genome-wide
# mapping of a quantitative trait. The data set is simulated assuming
# that all the genetic markers are uncorrelated with each other (i.e.,
# they are "unlinked").
#
# NOTE: Move this file later to the "demo" directory.
#
source("misc.R")
source("varbvs.R")
source("varbvsnorm.R")
source("varbvsnormupdate.R")
dyn.load("../src/diagsqr.so")
dyn.load("../src/varbvsr.so")

# SCRIPT PARAMETERS
# -----------------
n  <- 800   # Number of samples.
p  <- 2000  # Number of variables (genetic markers).
na <- 20    # Number of quantitative trait loci (QTLs).
se <- 4     # Variance of residual
r  <- 0.5   # Proportion of variance in trait explained by QTLs.

# Names of covariates.
covariates <- c("age","weight","glucose")

# Candidate values for the prior log-odds of inclusion.
logodds <- seq(-3,-1,0.1)

# Set the random number generator seed.
set.seed(1)

# GENERATE DATA SET
# -----------------
# Generate the minor allele frequencies so that they are uniform over range
# [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
# (i.e., unlinked), according to the specified minor allele frequencies.
cat("1. GENERATING DATA SET.\n")
maf <- 0.05 + 0.45*runif(p)
X   <- (runif(n*p) < maf) +
       (runif(n*p) < maf)
X   <- matrix(X,n,p,byrow = TRUE)
storage.mode(X) <- "double"

# Generate additive effects for the markers so that exactly na of them have
# a nonzero effect on the trait.
i       <- sample(p,na)
beta    <- rep(0,p)
beta[i] <- rnorm(na)

# Generate random labels for the markers.
colnames(X) <- paste0("rs",sample(1e6,p))

# Adjust the QTL effects so that we control for the proportion of variance
# explained (r). That is, we adjust beta so that r = a/(a+1), where I've
# defined a = beta'*cov(X)*beta. Here, sb is the variance of the (nonzero)
# QTL effects.
sb   <- r/(1-r)/var(c(X %*% beta))
beta <- sqrt(sb*se) * beta

# Generate a random intercept.
mu <- rnorm(1)

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

# Generate the quantitative trait measurements.
y <- mu + X %*% beta + sqrt(se)*rnorm(n)
if (m > 0)
  y <- y + Z %*% u
y <- c(y)

y <- y - c(Z %*% solve(crossprod(Z),c(y %*% Z)))
X <- X - Z %*% solve(crossprod(Z),t(Z) %*% X)
out <- varbvsnorm(X,y,1,1,rep(log(na/p),p),runif(p),rnorm(p))

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of a
# continuous outcome (quantitiative trait), with spike and slab priors on
# the coefficients.
cat("2. FITTING MODEL TO DATA.\n")
# TO DO.
# fit <- varbvs(X,Z,y,"gaussian",logodds = logodds)
