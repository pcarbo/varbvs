library(varbvs)
source("varbvssparse.R")

# SCRIPT PARAMETERS
# -----------------
k       <- 3           # Number of nonzero effects in sparse approx.
sigma   <- 9           # Variance of residual.
sa      <- 2           # Prior variance of nonzero effects.
logodds <- log(3/1000) # Prior inclusion probability.

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
cat("        variational    max.   incl variance params\n")
cat(" iter   lower bound  change   vars   sigma      sa\n")
alpha0  <- runif(p)
alpha0  <- alpha0/sum(alpha0)
mu0     <- rnorm(p)
logodds <- rep(logodds,1,p)
fit1    <- varbvsnorm(X,y,sigma,sa,logodds,alpha0,mu0,update.order = 1:p,
                      update.sigma = FALSE,update.sa = FALSE,tol = 1e-6)
cat("\n")

# FIT "SPARSE" VARIATIONAL APPROXIMATION
# --------------------------------------
cat("Fitting sparse variational approximation.\n")
fit2 <- varbvssparse(X,y,k,sigma,sa,logodds,tol = 1e-6)

