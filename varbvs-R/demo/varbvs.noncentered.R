# This script compares the "zero-centered" model, in which the prior
# on the regression coefficients is centered at zero, against the
# "non-centered" model, in which the prior mean is fitted to the data.
# The data set is simulated genotypes at unlinked genetic markers,
# with a continuous (i.e., quantitative) outcome.
library(lattice)
library(varbvs)  

# SCRIPT PARAMETERS
# -----------------
n  <- 800   # Number of samples.
p  <- 2000  # Number of variables (genetic markers).
na <- 20    # Number of quantitative trait loci (QTLs).
se <- 4     # Variance of residual.
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
X   <- matrix(as.double(X),n,p,byrow = TRUE)

# Generate non-negative additive effects for the markers so that
# exactly na of them have a nonzero effect on the trait.
i       <- sample(p,na)
beta    <- rep(0,p)
beta[i] <- abs(rnorm(na))

# Generate labels for the markers.
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

# Generate labels for the samples.
names(y)    <- sprintf("A%05d",sample(99999,n))
rownames(X) <- names(y)
if (!is.null(Z))
  rownames(Z) <- names(y)

# FIT ZERO-CENTERED MODEL
# -----------------------
cat("2. FITTING ZERO-CENTERED MODEL.\n")
fit1 <- varbvs(X,Z,y,"gaussian",logodds = logodds,n0 = 0,verbose = FALSE)

# FIT NON-CENTERED MODEL
# ----------------------
cat("3. FITTING NON-CENTERED MODEL.\n")
fit2 <- varbvs(X,Z,y,"gaussian",logodds = logodds,update.b0 = TRUE,
               n0 = 0,nb0 = 0,verbose = FALSE)

# COMPARE MODEL FIT
# -----------------
cat("Improvement in fit with non-centered model:\n")
cat(sprintf("Bayes factor = %0.2e\n",varbvsbf(fit1,fit2)))

# COMPARE ESTIMATES
# -----------------
# Plot the coefficients estimated using the centered and non-centered
# models. It is expected that the centered model will "shrink" the
# larger coefficients slightly more toward zero.
cat("4. PLOTTING COEFFICIENT ESTIMATES.\n")
trellis.par.set(par.xlab.text = list(cex = 0.75),
                par.ylab.text = list(cex = 0.75),
                axis.text = list(cex = 0.75))
markers <- labels(fit1)
b1      <- coef(fit1)
b1      <- b1[markers,ncol(b1)]
b2      <- coef(fit2)
b2      <- b2[markers,ncol(b2)]
print(xyplot(b2 ~ b1,data.frame(b1 = b1,b2 = b2),pch = 4,col = "black",
             cex = 0.6,
             panel = function(x, y, ...) {
               panel.xyplot(x,y,...)
               panel.abline(a = 0,b = 1,col = "magenta",lty = "dotted")
             },
             scales = list(limits = c(-1.1,1.1)),
             xlab = "coef (centered model)",
             ylab = "coef (non-centered model)"))
