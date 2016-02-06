# This script fits two variable selection models: the first ("null") model
# has a uniform prior for all variables (the 442,001 genetic markers), and
# the second model has higher prior probability for genetic markers near
# cytokine signaling genes. This analysis assess support for enrichment of
# Crohn's disease risk factors near cytokine signaling genes; a large Bayes
# factor means greater support for the enrichment hypothesis. The data in
# this analysis consist of 442,001 SNPs genotyped for 1,748 cases and 2,938
# controls.
library(varbvs)

# Initialize the random number generator. 
set.seed(1)

# LOAD GENOTYPES, PHENOTYPES AND PATHWAY ANNOTATION
# -------------------------------------------------
cat("LOADING DATA.\n")
load("/tmp/pcarbo/cd.RData")
load("~/data/cytokine.RData")

# FIT VARIATIONAL APPROXIMATION
# -----------------------------
# Compute the variational approximation given the assumption that all
# variables (genetic markers) are, *a priori*, equally likely to be included
# in the model.
cat("FITTING NULL MODEL TO DATA.\n")
fit_null <- varbvs(X,NULL,y,"binomial",logodds = -4)
