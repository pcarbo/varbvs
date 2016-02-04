# This script fits the Bayesian variable selection model to identify
# genetic markers associated with celiac disease risk. After removing
# all samples from the Finnish cohort, the data consist of 509,314
# SNPs genotyped for 3,149 cases and 6,325 controls, and principal
# components (PCs) which are used as covariates in the logistic
# regression.
library(varbvs)

# Initialize the random number generator. 
set.seed(1)

# LOAD GENOTYPE AND PHENOTYPE DATA
# --------------------------------
cat("LOADING DATA.\n")
load("/tmp/pcarbo/celiac.RData")

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a logistic regression model of
# a binary outcome (case-control status), with spike and slab priors
# on the coefficients.
cat("FITTING MODEL TO DATA.\n")
Z   <- as.matrix(panel[c("PC1","PC2")])
y   <- panel$pheno
fit <- varbvs(X,Z,y,"binomial",logodds = -5.5:0.25:-3)

# Compute "single-marker" posterior inclusion probabilities.
w   <- c(normalizelogweights(fit$logw))
pip <- varbvsindep(fit,X,Z,y) %*% w

# SAVE RESULTS
# ------------
cat("SAVING RESULTS.\n")
# TO DO.

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("SUMMARIZING RESULTS.\n")
varbvsprint(fit,n = 14)

# Show two "genome-wide scans", one using the posterior inclusion
# probabilities (PIPs) computed in the joint analysis of all
# variables, and one using the PIPs that ignore correlations between
# the variables. The latter is meant to look like a typical
# genome-wide "Manhattan" plot used to summarize the results of a
# genome-wide association study. Variables with PIP > 0.5 are
# highlighted.
i <- which(fit$alpha %*% w > 0.5)
