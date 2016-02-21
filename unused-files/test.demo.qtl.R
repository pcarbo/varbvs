library(varbvs)
library(glmnet)
source("varbvs-R/R/varbvs.R")
source("varbvs-R/R/varbvsnorm.R")
source("varbvs-R/R/misc.R")
source("varbvs-R/R/predict.varbvs.R")

# SCRIPT PARAMETERS
# -----------------
n  <- 1000  # Number of samples.
se <- 4     # Variance of residual.

# Set the random number generator seed.
set.seed(1)

# GENERATE DATA SET
# -----------------
# Generate the minor allele frequencies so that they are uniform over range
# [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
# (i.e., unlinked), according to the specified minor allele frequencies.
maf <- c(0.1,0.2)
X   <- (runif(n*2) < maf) +
       (runif(n*2) < maf)
X   <- matrix(as.double(X),n,2,byrow = TRUE)

# Specify additive effects.
beta <- c(-1,1)

# Specify the intercept.
mu <- 0.5

# Generate the covariate data (Z), and the linear effects of the
# covariates (u).
Z <- randn(n,3)
u <- c(-1,2,-3)

# Generate the quantitative trait measurements.
y <- c(mu + Z %*% u + X %*% beta + sqrt(se)*rnorm(n))

# Fit the elastic net model to the data.
fit.glm <- glmnet(cbind(Z,X),y,family = "gaussian",lambda = 0)

# Compute the estimates of Y using glmnet.
y.glm <- c(predict(fit.glm,cbind(Z,X)))

# Fit the model to the data.
fit <- varbvs(X,Z,y,"gaussian",logodds = seq(-1,0,0.1))

# Compute the estimates of Y using the variational approximation.
y.fit <- predict.varbvs(fit,X,Z)
