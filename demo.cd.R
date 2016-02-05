# This script fits the Bayesian variable selection model to identify genetic
# markers associated with Crohn's disease risk. The data consist of 442,001
# SNPs genotyped for 1,748 cases and 2,938 controls.
library(varbvs)
library(lattice)
library(latticeExtra)

# Initialize the random number generator. 
set.seed(1)

# LOAD GENOTYPE AND PHENOTYPE DATA
# --------------------------------
cat("LOADING DATA.\n")
load("/tmp/pcarbo/cd.RData")

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a logistic regression model of
# a binary outcome (case-control status), with spike and slab priors
# on the coefficients.
cat("FITTING MODEL TO DATA.\n")
fit <- varbvs(X,NULL,y,"binomial",logodds = seq(-6,-3,0.25))

# Compute "single-marker" posterior inclusion probabilities.
w   <- c(normalizelogweights(fit$logw))
pip <- varbvsindep(fit,X,NULL,y) %*% w

# SAVE RESULTS
# ------------
cat("SAVING RESULTS.\n")
save(list = c("fit","map","pip"),file = "/tmp/pcarbo/varbvs.demo.cd.RData")

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("SUMMARIZING RESULTS.\n")
varbvsprint(fit,n = 9)

# Show two "genome-wide scans", one using the posterior inclusion
# probabilities (PIPs) computed in the joint analysis of all
# variables, and one using the PIPs that ignore correlations between
# the variables. The latter is meant to look like a typical
# genome-wide "Manhattan" plot used to summarize the results of a
# genome-wide association study. Variables with PIP > 0.5 are
# highlighted.
trellis.device(height = 3,width = 10)
i <- which(fit$alpha %*% w > 0.5)
varbvsplot(fit,groups = map$chr,vars = i,gap = 7500,xlab = "posterior prob.")

