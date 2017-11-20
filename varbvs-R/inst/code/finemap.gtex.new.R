library(varbvs)
source("varbvsfinemap.R")

# SCRIPT PARAMETERS
# -----------------
sigma   <- 1            # Variance of residual.
sa      <- 1            # Prior variance of nonzero effects.
logodds <- log(3/1000)  # Prior inclusion probability.

# These are the additive effects used to simulate the
# phenotype data.
beta                 <- rep(0,1000)
beta[c(201,562,952)] <- c(1,-1,1)

# LOAD GENOTYPE DATA
# ------------------
# Retrieve genotype data for 1,000 SNPs near gene FMO2.
cat("Loading genotype data.\n")
load("../datafiles/Thyroid.FMO2.1Mb.RData")
rm(Z,y)
X <- X[,3501:4500]
n <- nrow(X)
p <- ncol(X)
storage.mode(X) <- "double"

# Center columns of X.
X <- scale(X,center = TRUE,scale = FALSE)

# SIMULATE PHENOTYPE DATA
# -----------------------
cat("Simulating phenotype data.\n")
y <- c(X %*% beta + sqrt(sigma)*rnorm(n))
y <- y - mean(y)

# FIT FULLY-FACTORIZED VARIATIONAL APPROXIMATION
# ----------------------------------------------
cat("Fitting fully-factorized variational approximation.\n")
alpha0 <- runif(p)
alpha0 <- alpha0/sum(alpha0)
mu0    <- rnorm(p)
fit1   <- varbvsnorm(X,y,sigma,sa,logodds,alpha0,mu0,
                     update.sigma = FALSE,update.sa = FALSE)

