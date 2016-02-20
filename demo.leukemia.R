# TO DO: Explain here what this script does, and how to use it. Also
# cite Golub et al, elastic net paper, and Friedman et al (2010)
# J. Stat. Soft. paper.
library(glmnet)
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
# TO DO.

# LOAD LEUKEMIA DATA
# ------------------
cat("1. Loading leukemia data.\n")
load("varbvs-R/data/leukemia.RData")
X <- Leukemia$x
y <- Leukemia$y

# Set the random number generator seed.
set.seed(1)

# FIT ELASTIC NET MODEL TO DATA
# -----------------------------
cat("2. Running 20-fold cross-validation to select L1-penalty strength.\n")
out <- cv.glmnet(X,y,family = "binomial",alpha = 0.2,nfolds = 20,
                 type.measure = "class")

# Show the misclassification error at different settings of the
# L1-penalty strength. Then choose the largest value of the
# L1-penalty strength (lambda) that is within 1 standard error of the
# smallest misclassification error.
dev.new(height = 3,width = 7)
plot(out)
lambda <- out$lambda.1se

# Fit the elastic net model to the data.
cat(sprintf("3. Fitting elastic net model with lambda = %0.3f.\n",lambda))
fit.glmnet <- glmnet(X,y,family = "binomial",alpha = 0.2)
dev.new(height = 3,width = 5)
plot(fit.glmnet,xvar = "lambda",xlim = c(-2,1))

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a logistic regression model of
# the binary outcome (type of leukemia), with spike-and-slab priors
# on the coefficients.
cat("4. Fitting Bayesian variable selection model to data.\n")
fit.varbvs <- varbvs(X,NULL,y,"binomial",logodds = seq(-4,-1,0.1),sa = 0.2)

