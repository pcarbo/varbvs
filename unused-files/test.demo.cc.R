library(glmnet)
library(varbvs)  
source("varbvs-R/R/varbvs.R")
source("varbvs-R/R/varbvsbin.R")
source("varbvs-R/R/varbvsbinz.R")
source("varbvs-R/R/misc.R")
source("varbvs-R/R/predict.varbvs.R")

# SCRIPT PARAMETERS
# -----------------
n <- 5000

# Set the random number generator seed.
set.seed(1)

# GENERATE DATA SET
# -----------------
# Generate the minor allele frequencies so that they are uniform over range
# [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
# (i.e., unlinked), according to the specified minor allele frequencies.
maf <- c(0.25,0.25)
X   <- (runif(n*2) < maf) +
       (runif(n*2) < maf)
X   <- matrix(as.double(X),n,2,byrow = TRUE)

# Specify additive effects.
beta <- c(-1,1)

# Specify the intercept.
mu <- 0.5

# Generate the covariate data (Z), and the linear effects of the
# covariates (u).
Z <- randn(n,2)
u <- c(-2,2)

# For each sample, calculate the probability of being a case (y = 1).
w <- c(mu + Z %*% u + X %*% beta)

# Simulate the binary trait (case-control status) as a coin toss with
# success rates given by the logistic regression.
y <- as.double(runif(n) < sigmoid(w))

# Fit the elastic net model to the data.
fit.glm <- glmnet(cbind(Z,X),y,family = "binomial",lambda = 0)

# Compute the estimates of Y using glmnet.
y.glm <- as.integer(c(predict(fit.glm,cbind(Z,X),type = "class")))

# Fit the variable selection model to the data.
fit <- varbvs(X,Z,y,"binomial",logodds = seq(-1,0,0.1))

# Compute the estimates of Y using varbvs.
y.fit <- predict.varbvs(fit,X,Z)

# Summarize glmnet predictions.
print(table(factor(y),factor(y.glm)))

# Summarize varbvs predictions.
print(table(factor(y),factor(y.fit)))
