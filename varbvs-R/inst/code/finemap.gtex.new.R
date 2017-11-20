library(varbvs)
source("varbvsfinemap.R")

# SCRIPT PARAMETERS
# -----------------
# TO DO.

# LOAD GENOTYPE DATA
# ------------------
# Retrieve genotype data for 1,000 SNPs near gene FMO2.
cat("Loading genotype data.\n")
load("../datafiles/Thyroid.FMO2.1Mb.RData")
storage.mode(X) <- "double"
X <- X[,3501:4500]
n <- nrow(X)
p <- ncol(X)

# SIMULATE PHENOTYPE DATA
cat("Simulating phenotype data.\n")

varbvs:::varbvsnormupdate(X,sigma,sa, logodds, xy, d, alpha0, mu0, Xr0, i,
              algorithm.version = c(".Call","Rcpp","R")) {
