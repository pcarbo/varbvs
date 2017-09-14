# TO DO: Explain here what this script does, and how to use it.
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
n   <- 800   # Number of samples.
p   <- 2000  # Number of variables (SNPs).
r   <- 0.5   # Proportion of variance in trait explained by SNPs.
rho <- 0.3   # Proportion of additive genetic variance due to large effects.
na  <- 20    # Number of large (QTL) effects.

# Candidate values for the prior log-odds of inclusion.
logodds <- seq(-3,-1,0.1)

# Candidate values for the log-odds of being a QTL (logodds), the
# proportion of variance explained (h), and the proportion of additive
# genetic variance due to "large" QTL effects (d).
logodds <- seq(-3,-1,0.25)
h       <- seq(0.1,0.9,0.1)
d       <- seq(0.1,0.9,0.1)

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

# Generate the "small" (polygenic) effects for the SNPs.
u <- randn(p,1)

# Generate the "large" (QTL) additive effects for the SNPs.
i    <- sample(p,na)
b    <- rep(0,p)
b[i] <- rnorm(na)

# Generate random labels for the markers.
colnames(X) <- paste0("rs",sample(1e6,p))

# Adjust the additive effects so that we control for (1) the
# proportion of additive genetic variance that is due to QTL effects
# (rho), and (2) the proportion of variance explained (r). That is, we
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
st   <- r/(1-r) * d/var(c(X %*% b))
beta <- sqrt(st) * b
if (d == 0)
   case 0
    sa = r/(1-r)/var(X*u,1);
   case 1
    sa = 0;
   otherwise
    sa = max(roots2(var(X*u,1),2*dot(X*beta,X*u)/n,var(X*beta,1) - r/(1-r)))^2;
  end
  u = sqrt(sa) * u;
    
  % Generate the quantitative trait measurements.
  y = X*(u + beta) + randn(n,1);
