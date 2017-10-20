# This is similar to demo.mix.R, except that the prior is a larger
# mixture of normals. This script is mainly intended to illustrate the
# use of the "drop.threshold" varbvsmix argument to speed up
# computation when a lot of the mixture components have a negligible
# (near zero) weight.
library(lattice)
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
sd <- c(0,   0.1, 0.2, 0.5)
w  <- c(0.95,0.03,0.01,0.01)

# The dense grid of prior standard deviations used to model the data.
sd.grid <- c(0,10^seq(-2,1,length.out = 19))

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
fit <- varbvsmix(X,Z,y,sd.grid^2)

# Plot the estimated coefficients against the ground-truth coefficients.
trellis.par.set(par.xlab.text = list(cex = 0.75),
                par.ylab.text = list(cex = 0.75),
                axis.text = list(cex = 0.75))
print(xyplot(beta.est ~ beta.true,
             data.frame(beta.true = beta,
                        beta.est  = rowSums(fit$alpha * fit$mu)),
             pch = 4,col = "black",cex = 0.6,
             panel = function(x, y, ...) {
               panel.xyplot(x,y,...)
               panel.abline(a = 0,b = 1,col = "magenta",lty = "dotted")
             },
             scales = list(limits = c(-1.1,1.1)),
             xlab = "ground-truth regression coefficient",
             ylab = "estimated regression coefficient"),
      split = c(1,1,3,1),
      more = TRUE)

# Show the change in the variational lower bound at each iteration of the
# co-ordinate ascent algorithm.
numiter <- length(fit$logZ)
print(xyplot(y ~ x,data.frame(x = 1:numiter,y = max(fit$logZ) - fit$logZ),
             type = "l",col = "darkorange",lwd = 2,
             scales = list(y = list(log = 10)),xlab = "iteration",
             ylab = "distance from final lower bound"),
      split = c(2,1,3,1),
      more = TRUE)

# Plot the number of nonzero mixture components at each iteration of
# the co-ordinate ascent algorithm.
print(xyplot(y ~ x,data.frame(x = 1:numiter,y = length(sd.grid) - fit$nzw),
             type = "l",col = "royalblue",lwd = 2,xlab = "iteration",
             ylab = "nonzero mixture components"),
      split = c(3,1,3,1),
      more = TRUE)
