# This script fits the Bayesian variable selection model to identify
# genetic markers associated with Crohn's disease risk. The data
# consist of 442,001 SNPs genotyped for 1,748 cases and 2,938 controls.
#
# Note that file cd.RData cannot be made publicly available due to
# data sharing restrictions, so this script is for viewing only.
#
library(varbvs)
library(lattice)

# Initialize the random number generator. 
set.seed(1)

# LOAD GENOTYPE AND PHENOTYPE DATA
# --------------------------------
cat("LOADING DATA.\n")
load("cd.RData")

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a logistic regression model of
# a binary outcome (case-control status), with spike and slab priors
# on the coefficients.
cat("FITTING MODEL TO DATA.\n")
r <- system.time(fit <- varbvs(X,NULL,y,"binomial",logodds = seq(-6,-3,0.25)))
cat(sprintf("Modeling fitting took %0.2f minutes.\n",r["elapsed"]/60))

# Compute "single-marker" posterior inclusion probabilities.
w   <- c(normalizelogweights(fit$logw))
pip <- c(varbvsindep(fit,X,NULL,y)$alpha %*% w)

# SAVE RESULTS
# ------------
cat("SAVING RESULTS.\n")
save(list = c("fit","map","pip","r"),file = "varbvs.demo.cd.RData")

# SUMMARIZE POSTERIOR DISTRIBUTION
# --------------------------------
cat("SUMMARIZING RESULTS.\n")
print(summary(fit,nv = 9))

# Show two "genome-wide scans", one using the posterior inclusion
# probabilities (PIPs) computed in the joint analysis of all
# variables, and one using the PIPs that ignore correlations between
# the variables. The latter is meant to look like a typical
# genome-wide "Manhattan" plot used to summarize the results of a
# genome-wide association study. Variables with PIP > 0.5 are
# highlighted.
trellis.device(height = 4,width = 10)
i <- which(fit$alpha %*% w > 0.5)
var.labels <- paste0(round(map$pos[i]/1e6,digits = 2),"Mb")
print(plot(fit,groups = map$chr,vars = i,var.labels = var.labels,gap = 7500,
           ylab = "posterior prob."),
      split = c(1,1,1,2),more = TRUE)
print(plot(fit,groups = map$chr,score = log10(pip + 0.001),vars = i,
           var.labels = var.labels,gap = 7500,ylab = "log10 posterior prob."),
      split = c(1,2,1,2),more = FALSE)
