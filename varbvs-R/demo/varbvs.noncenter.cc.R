# This script compares the "zero-centered" model, in which the prior
# on the regression coefficients is centered at zero, against the
# "non-centered" model, in which the prior mean is fitted to the
# data. The data set is simulated genotypes at unlinked genetic
# markers, with a binary (e.g., case-control) outcome.
library(lattice)
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
n  <- 1400  # Number of samples (subjects).
p  <- 1000  # Number of variables (genetic markers).
na <- 20    # Number of markers that affect the binary outcome.
sa <- 0.25  # Variance of log-odds ratios.
p1 <- 0.01  # Proportion of subjects that will be cases when all the
            # genotypes are zero.

# Names of covariates.
covariates <- c("age","weight")

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

# Generate non-negative additive effects for the markers so that
# exactly na of them have a positive effect on the trait.
i       <- sample(p,na)
beta    <- rep(0,p)
beta[i] <- sqrt(sa)*abs(rnorm(na))

# Generate labels for the markers.
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
mu <- varbvs:::logit(p1)
w  <- mu + X %*% beta
if (m > 0)
  w <- w + Z %*% u

# Simulate the binary trait (case-control status) as a coin toss with
# success rates given by the logistic regression.
y <- as.double(runif(n) < varbvs:::sigmoid(w))

# Generate labels for the samples.
names(y)    <- sprintf("A%05d",sample(99999,n))
rownames(X) <- names(y)
if (!is.null(Z))
  rownames(Z) <- names(y)

# FIT ZERO-CENTERED MODEL
# -----------------------
cat("2. FITTING ZERO-CENTERED MODEL.\n")
fit1 <- varbvs(X,Z,y,"binomial",logodds = logodds,n0 = 0,verbose = TRUE)

# FIT NON-CENTERED MODEL
# ----------------------
cat("3. FITTING NON-CENTERED MODEL.\n")
fit2 <- varbvs(X,Z,y,"binomial",logodds = logodds,update.sa = FALSE,
               update.b0 = TRUE,n0 = 0,nb0 = 0,verbose = TRUE)

# COMPARE MODEL FIT
# -----------------
cat("Improvement in fit with non-centered model:\n")
cat(sprintf("Bayes factor = %0.2e\n",varbvsbf(fit1,fit2)))

# EVALUATE MODEL PREDICTIONS
# --------------------------
# Evaluate accuracy of the two fitted models.
cat("4. EVALUATING FITTED MODELS.\n")
y1 <- predict(fit1,X,Z,type = "class")
y2 <- predict(fit2,X,Z,type = "class")
cat("r^2 between predicted Y and observed Y\n")
cat(sprintf("zero-centered model: %0.3f\n",cor(y,y1)^2))
cat(sprintf("non-centered model:  %0.3f\n",cor(y,y2)^2))

# COMPARE ESTIMATES
# -----------------
# Plot the coefficients estimated using the centered and non-centered
# models. It is expected that the centered model will "shrink" the
# larger coefficients slightly more toward zero.
cat("5. PLOTTING COEFFICIENT ESTIMATES.\n")
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
