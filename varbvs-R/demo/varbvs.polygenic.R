# TO DO: Explain here what this script does.
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
n  <- 800   # Number of samples.
p  <- 2000  # Number of variables (genetic markers).
na <- 20    # Number of quantitative trait loci (QTLs).
se <- 4     # Variance of residual.
r  <- 0.5   # Proportion of variance in trait explained by markers.
d  <- 0.6   # Proportion of additive genetic variance due to QTLs.

# Names of covariates.
covariates <- c("age","weight")

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
X   <- matrix(as.double(X),n,p,byrow = TRUE)

# Generate (small) polygenic additive effects for the markers.
u <- rnorm(p)

# Generate (large) QTL effects for the markers.
i     <- sample(p,na)
bx    <- rep(0,p)
bx[i] <- rnorm(na)

# Generate labels for the markers.
colnames(X) <- paste0("rs",sample(1e6,p))

# Adjust the additive effects so that we control for the proportion of
# additive genetic variance that is due to QTL effects (d) and the
# total proportion of variance explained (r). That is, we adjust beta
# and u so that
#
#   r = a/(a+1)
#   d = b/a,
#
# where I've defined
#
#   a = (u + beta)'*cov(X)*(u + beta),
#   b = beta'*cov(X)*beta.
#
# Note: this code only works if d or r are not exactly 0 or exactly 1.
st <- c(r/(1-r) * d/var(X %*% bx))
bx <- sqrt(st) * bx
sa <- max(Re(polyroot(c(c(var(X %*% bx) - r/(1-r)),
                        2*sum((X %*% bx) * (X %*% u))/n,
                        c(var(X %*% u))))))^2
u    <- sqrt(sa) * u

# Generate a random intercept.
mu <- rnorm(1)

# Generate the covariate data (Z), and the linear effects of the
# covariates (u).
m <- length(covariates)
if (m > 0) {
  Z  <- randn(n,m)
  bz <- rnorm(m)
  colnames(Z) <- covariates
} else {
  Z <- NULL
}

# Generate the quantitative trait measurements.
y <- mu + X %*% (u + bx) + sqrt(se)*rnorm(n)
if (m > 0)
  y <- y + Z %*% u
y <- c(y)

# Generate labels for the samples.
names(y)    <- sprintf("A%05d",sample(99999,n))
rownames(X) <- names(y)
if (!is.null(Z))
  rownames(Z) <- names(y)
