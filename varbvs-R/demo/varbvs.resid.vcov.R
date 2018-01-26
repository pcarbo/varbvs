# TO DO: Explain here what this script does.
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
n  <- 800   # Number of samples.
p  <- 2000  # Number of variables (genetic markers).
na <- 2     # Number of quantitative trait loci (QTLs).
r  <- 0.5   # Proportion of variance in trait explained by markers.
d  <- 0.6   # Proportion of additive genetic variance due to QTLs.

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
X   <- scale(X,center = TRUE,scale = FALSE)

# Generate (small) polygenic additive effects for the markers.
u <- rnorm(p)

# Generate (large) QTL effects for the markers.
i       <- sample(p,na)
beta    <- rep(0,p)
beta[i] <- rnorm(na)

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
st   <- c(r/(1-r) * d/var(X %*% beta))
beta <- sqrt(st) * beta
sa   <- max(Re(polyroot(c(c(var(X %*% beta) - r/(1-r)),
                          2*sum((X %*% beta) * (X %*% u))/n,
                          c(var(X %*% u))))))^2
u    <- sqrt(sa) * u

# Generate the quantitative trait measurements.
y <- X %*% (u + beta) + rnorm(n)
y <- c(y)

# Generate labels for the samples.
names(y)    <- sprintf("A%05d",sample(99999,n))
rownames(X) <- names(y)

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------

# TO DO: Explain what this is doing.
cat("2. FITTING MODEL TO DATA.\n")
fit1 <- varbvs(X,NULL,y,"gaussian",logodds = logodds,n0 = 0)
fit2 <- varbvs(X,NULL,y,"gaussian",logodds = logodds,n0 = 0,
               resid.vcov = diag(n) + sa*tcrossprod(X))
cat("\n")


# SUMMARIZE RESULTS
# -----------------
cat("3. SUMMARIZING RESULTS.\n")
print(summary(fit1))
cat("\n")
