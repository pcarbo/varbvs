# The varbvs and varbvsmix functions should produce the same estimates
# when there are exactly 2 mixture components (a "spike" and a
# "slab"). This script verifies this in a small simulated data set.
library(varbvs)  

# SCRIPT PARAMETERS
# -----------------
n  <- 1000  # Number of samples.
p  <- 2000  # Number of variables (genetic markers).
se <- 4     # Variance of residual.

# Names of the covariates.
covariates <- c("age","weight")

# The standard deviations and mixture weights used to simulate the
# additive effects on the quantitative trait. Note that the first
# mixture component must have a standard deviation of exactly zero.
sd <- c(0,0.5)
w  <- c(0.95,0.05)

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
# Fit the varbvsmix model.
cat("2. FITTING VARBVSMIX MODEL TO DATA.\n")
alpha   <- runif(p)
alpha   <- cbind(alpha,1 - alpha)
mu      <- cbind(0,rnorm(p))
fit     <- varbvsmix(X,Z,y,sd^2,alpha = alpha,mu = mu)

# Fit the varbvs model.
cat("3. FITTING VARBVS MODEL TO DATA.\n")
w1   <- fit$w[2]
fit2 <- varbvs(X,Z,y,sa = sd[2]^2,logodds = log10(w1/(1 - w1)),
               alpha = alpha[,2],mu = mu[,2])

# Check that the parameter estimates are the same.
niter <- length(fit$logZ)
relerr <- function (x, y)
  abs(x - y)/min(abs(x),abs(y))
cat(sprintf("Max. difference in PIPs:             %0.2e\n",
            max(abs(fit$alpha[,2] - fit2$alpha))))
cat(sprintf("Max. difference in posterior means:  %0.2e\n",
            max(abs(fit$mu[,2] - fit2$mu))))
cat(sprintf("Relative difference in lower bounds: %0.2e\n",
            relerr(fit$logZ[niter],fit2$logw)))
