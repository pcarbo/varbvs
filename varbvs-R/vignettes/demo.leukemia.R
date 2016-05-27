# This script demonstrates application of *glmnet* and *varbvs* to the
# Leukemia data set. The main aim of this script is to illustrate some
# of the different properties of Bayesian variable selection and
# penalized sparse regression (as implemented by varbvs and glmnet,
# respectively). This script also reproduces the results and graphs
# presented in the first example of Carbonetto et al (2016).
#
# We use the preprocessed data of Dettling (2004) retrieved from the
# supplementary materials accompanying Friedman et al (2010). The data
# are represented as a 72 x 3,571 matrix of gene expression values
# (variable X), and a vector of 72 binary disease outcomes (variable
# y).
#
# REFERENCES
#
# Carbonetto, P., Zhou, X., Trynka, G., Stephens, M. (2016) Fast
# variable selection for genome-wide association studies and other
# large-scale regression applications. *Forthcoming*.
#
# Dettling, M. (2004). BagBoosting for tumor classification with gene
# expression data. Bioinformatics 20, 3583–3593.
#
# Friedman, J., Hastie, T., Tibshirani, R. (2010) Regularization paths
# for generalized linear models via coordinate descent. Journal of
# Statistical Software 33, 1–22.
#
library(lattice)
library(latticeExtra)
library(glmnet)
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
# glmnet settings.
nfolds <- 20                  # Number of cross-validation folds.
alpha  <- 0.95                # Elastic net mixing parameter.
lambda <- 10^(seq(-2,0,0.05)) # Lambda sequence.

# varbvs settings.
logodds <- seq(-3.5,-1.5,0.1) # Candidate prior log-odds settings.
sa      <- 1                  # Prior variance of coefficients.

# LOAD LEUKEMIA DATA
# ------------------
cat("1. Loading leukemia data.\n")
data(leukemia)
X <- leukemia$x
y <- leukemia$y
rm(leukemia)

# Set the random number generator seed.
set.seed(1)

# FIT ELASTIC NET MODEL TO DATA
# -----------------------------
# First run 20-fold cross-validation to select the largest setting of
# lambda that is within one standard error of the minimum
# classification error.
cat("2. Running 20-fold cross-validation to select L1-penalty strength.\n")
r <- system.time(out.cv.glmnet <-
       cv.glmnet(X,y,family = "binomial",type.measure = "class",
                 alpha = alpha,nfolds = nfolds,lambda = lambda))
lambda <- out.cv.glmnet$lambda
cat(sprintf("Cross-validation took %0.2f seconds.\n",r["elapsed"]))
rm(r)

# Show the classification error at different settings of lambda. Then
# choose the largest value of lambda that is within 1 standard error
# of the smallest misclassification error.
lambda.opt <- out.cv.glmnet$lambda.1se

# Fit the elastic net model to the data.
cat("3. Fitting elastic net model.\n")
r <- system.time(fit.glmnet <-
       glmnet(X,y,family = "binomial",lambda = lambda,alpha = alpha))
cat(sprintf("Model fitting took %0.2f seconds.\n",r["elapsed"]))
rm(r)
                 
# Compute estimates of the disease outcome using the fitted model, and
# compare against the observed values.
cat("4. Summarizing results of glmnet analysis.\n")
cat("classification results with lambda = ",lambda.opt,":\n",sep="")
y.glmnet <- as.integer(c(predict(fit.glmnet,X,s = lambda.opt,type = "class")))
print(table(true = factor(y),pred = factor(y.glmnet)))

# Plot the classification error at different settings of lambda.
trellis.device(height = 4.5,width = 7)
trellis.par.set(par.xlab.text = list(cex = 0.65),
                par.ylab.text = list(cex = 0.65),
                axis.text     = list(cex = 0.65))
Y       <- predict(fit.glmnet,X,type = "class")
mode(Y) <- "numeric"
print(with(out.cv.glmnet,
           xyplot(y ~ x,data.frame(x = log10(lambda),y = cvm),type = "l",
                  col = "blue",xlab = "log10 lambda",
                  ylab = "classification error",
                  scales = list(y = list(limits = c(-0.02,0.45))),
                  panel = function(x, y, ...) {
                    panel.xyplot(x,y,...)
                    panel.abline(v = log10(lambda.opt),col = "orangered",
                                 lwd = 1,lty = "dotted")
                  }) +
           as.layer(xyplot(y ~ x,data.frame(x = log10(lambda),y = cvm),
                           pch = 20,cex = 0.6,col = "blue")) +
           as.layer(xyplot(y ~ x,data.frame(x = log10(lambda),y = cvup),
                           type = "l",col = "blue",lty = "dashed")) +
           as.layer(xyplot(y ~ x,data.frame(x = log10(lambda),y = cvlo),
                           type = "l",col = "blue",lty = "dashed")) +
           as.layer(xyplot(y ~ x,data.frame(x = log10(lambda),
                                            y = colMeans(abs(Y - y))),
                           type = "l",col = "red",lwd = 2))),
           split = c(1,1,2,2),more = TRUE)
rm(Y)

# Show the number of nonzero regression coefficients at different
# settings of lambda.
print(with(out.cv.glmnet,
           xyplot(y ~ x,data.frame(x = log10(lambda),y = nzero),type = "l",
                  col = "blue",xlab = "log10 lambda",
                  ylab = "num. nonzero coefs",
                  panel = function(x, y, ...) {
                    panel.xyplot(x,y,...)
                    panel.abline(v = log10(lambda.opt),col = "orangered",
                                 lwd = 1,lty = "dotted")
                  }) +
           as.layer(xyplot(y ~ x,data.frame(x = log10(lambda),y = nzero),
                           pch = 20,cex = 0.6,col = "blue"))),
      split = c(1,2,2,2),more = TRUE)

# Plot evolution of regression coefficients at different settings of
# lambda. Do not show the intercept.
vars <- setdiff(which(rowSums(abs(coef(fit.glmnet))) > 0),1)
n    <- length(vars)
b    <- as.matrix(t(coef(fit.glmnet)[vars,]))
r    <- xyplot(y ~ x,data.frame(x = log10(lambda),y = b[,1]),type = "l",
               col = "blue",xlab = "log10 lambda",
               ylab = "regression coefficient",
               scales = list(x = list(limits = c(-2.35,0.1)),
                             y = list(limits = c(-0.8,1.2))),
               panel = function(x, y, ...) {
                 panel.xyplot(x,y,...);
                 panel.abline(v = log10(lambda.opt),col = "orangered",
                              lwd = 1,lty = "dotted");
                 ltext(x = -2,y = b[nrow(b),],labels = colnames(b),pos = 2,
                       offset = 0.5,cex = 0.5)
               })
for (i in 2:n)
  r <- r + as.layer(xyplot(y ~ x,data.frame(x = log10(lambda),y = b[,i]),
                           type = "l",col = "blue"))
print(r,split = c(2,1,2,1),more = FALSE)
rm(vars,n,b,r,i)

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a logistic regression model of
# the binary outcome (type of leukemia), with spike-and-slab priors
# on the coefficients.
cat("5. Fitting Bayesian variable selection model to data.\n")
r <- system.time(fit.varbvs <-
       varbvs(X,NULL,y,"binomial",logodds = logodds,sa = sa,verbose = FALSE))
cat(sprintf("Model fitting took %0.2f seconds.\n",r["elapsed"]))
rm(r)

# Compute estimates of the disease outcome using the fitted model, and
# compare against the observed values.
cat("6. Summarizing results of varbvs analysis.\n")
cat("classification results from fitted varbvs model:\n")
y.varbvs <- predict(fit.varbvs,X)
print(table(true = factor(y),pred = factor(y.varbvs)))

# Show the classification error at each setting of the prior log-odds.
trellis.device(height = 4.5,width = 7)
trellis.par.set(par.xlab.text = list(cex = 0.65),
                par.ylab.text = list(cex = 0.65),
                axis.text     = list(cex = 0.65))
m    <- length(logodds)
err  <- rep(0,m)
logw <- fit.varbvs$logw
for (i in 1:m) {
  fit.varbvs$logw[]  <- 0
  fit.varbvs$logw[i] <- 100
  ypred  <- predict(fit.varbvs,X)
  err[i] <- mean(y != ypred)
}
fit.varbvs$logw <- logw
print(xyplot(y ~ x,data.frame(x = logodds,y = err),type = "l",
             col = "blue",xlab = "prior log-odds",
             ylab = "classification error",
             scales = list(x = list(limits = c(-1.4,-3.6)))) +
      as.layer(xyplot(y ~ x,data.frame(x = logodds,y = err),
                      col = "blue",pch = 20,cex = 0.65)),
      split = c(1,1,2,2),more = TRUE)
rm(err,logw,ypred,i)

# Show probability density of prior log-odds.
w        <- normalizelogweights(fit.varbvs$logw)
names(w) <- logodds
print(xyplot(y ~ x,data.frame(x = logodds,y = w),type = "l",col = "blue",
             xlab = "prior log-odds",ylab = "posterior prob.",
             scales = list(x = list(limits = c(-1.4,-3.6)))) +
      as.layer(xyplot(y ~ x,data.frame(x = logodds,y = w),
                      col = "blue",pch = 20,cex = 0.65)),
      split = c(1,2,2,2),more = TRUE)
rm(w)

# Plot evolution of posterior inclusion probabilities (PIPs) at
# different settings of the prior log-odds.
n     <- 10
vars  <- order(fit.varbvs$alpha[,m],decreasing = TRUE)[1:n]
alpha <- t(fit.varbvs$alpha[vars,] * fit.varbvs$mu[vars,])
alpha <- pmax(alpha,-0.075)
r     <- xyplot(y ~ x,data.frame(x = logodds,y = (alpha[,1])),
                type = "l",col = "blue",xlab = "prior log-odds",
                ylab = "regression coefficient",
                scales = list(x = list(limits = c(-1.15,-3.6)),
                              y = list(limits = c(-0.08,0.06))),
                panel = function(x, y, ...) {
                  panel.xyplot(x,y,...);
                  ltext(x = -1.5,y = (alpha[m,]),
                        labels = colnames(alpha),pos = 2,
                        offset = 0.5,cex = 0.5)
                })
for (i in 2:n)
  r <- r + as.layer(xyplot(y ~ x,
                           data.frame(x = logodds,y = (alpha[,i])),
                           type = "l",col = "blue"))
print(r,split = c(2,1,2,1),more = FALSE)
rm(m,n,i,vars,alpha,r)


