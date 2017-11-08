
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
na <- 10    # Number of markers that affect the binary outcome.
sa <- 0.2   # Variance of log-odds ratios.
p1 <- 0.1   # Proportion of subjects that will be cases when all the
            # genotypes are zero.

# Names of covariates.
# covariates <- c("age","weight")
covariates <- NULL

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

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the zero-centered model.
cat("2. FITTING ZERO-CENTERED MODEL.\n")
fit1 <- varbvs(X,Z,y,"binomial",b0 = 0,logodds = logodds,n0 = 0,
               verbose = FALSE)

# Fit the non-zero-centered model.
cat("3. FITTING NON-ZERO-CENTERED MODEL.\n")
fit2 <- varbvs(X,Z,y,"binomial",logodds = logodds,update.b0 = TRUE,
               n0 = 0,nb0 = 0,verbose = FALSE)

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("4. SUMMARIZING NON-ZERO-CENTERED MODEL.\n")
print(summary(fit2))
cat("\n")

# COMPARE ESTIMATES AGAINST GROUND-TRUTH
# --------------------------------------
# Plot the estimated coefficients against the ground-truth coefficients.
# It is expected that oefficients near zero will be "shrunk" to zero.
cat("5. PLOTTING COEFFICIENT ESTIMATES.\n")
trellis.par.set(par.xlab.text = list(cex = 0.75),
                par.ylab.text = list(cex = 0.75),
                axis.text = list(cex = 0.75))
coef.scatterplot <- function (fit, beta.true) {
  markers  <- labels(fit)
  beta.est <- coef(fit)
  beta.est <- beta.est[markers,ncol(beta.est)]
  return(xyplot(beta.est ~ beta.true,
                data.frame(beta.true = beta,beta.est = beta.est),
                pch = 4,col = "black",cex = 0.6,
                panel = function(x, y, ...) {
                  panel.xyplot(x,y,...)
                  panel.abline(a = 0,b = 1,col = "magenta",lty = "dotted")
                },
                scales = list(limits = c(-0.1,1.1)),
                xlab = "ground-truth regression coefficient",
                ylab = "estimated regression coefficient"))
}
print(coef.scatterplot(fit1,coef(fit2)),split = c(1,1,2,1),more = TRUE)
print(coef.scatterplot(fit2),split = c(2,1,2,1),more = FALSE)

# EVALUATE MODEL PREDICTIONS
# --------------------------
# Evaluate accuracy of the two fitted models.
cat("6. EVALUATING FITTED MODELS.\n")
y1 <- predict(fit1,X,Z,type = "class")
y2 <- predict(fit2,X,Z,type = "class")
cat("r^2 between predicted Y and observed Y\n")
cat(sprintf("zero-centered model: %0.3f\n",cor(y,y1)^2))
cat(sprintf("non-centered model:  %0.3f\n",cor(y,y2)^2))

