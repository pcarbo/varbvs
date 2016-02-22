# TO DO: Explain here what this script does, and how to use it. Also
# cite Golub et al, elastic net paper, and Friedman et al (2010)
# J. Stat. Soft. paper.
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
sa      <- 1                   # Prior variance of coefficients.
logodds <- seq(-3.5,-1.5,0.1)  # Candidate prior log-odds settings.

# LOAD LEUKEMIA DATA
# ------------------
cat("1. Loading leukemia data.\n")
load("varbvs-R/data/leukemia.RData")
X <- Leukemia$x
y <- Leukemia$y
rm(Leukemia)

# Set the random number generator seed.
set.seed(1)

# FIT ELASTIC NET MODEL TO DATA
# -----------------------------
# First run 20-fold cross-validation to select the largest setting of
# lambda that is within one standard error of the minimum
# classification error.
cat("2. Running 20-fold cross-validation to select L1-penalty strength.\n")
out.cv.glmnet <-
  cv.glmnet(X,y,family = "binomial",type.measure = "class",alpha = alpha,
            nfolds = nfolds,lambda = lambda)

# Show the classification error at different settings of lambda. Then
# choose the largest value of lambda that is within 1 standard error
# of the smallest misclassification error.
lambda.opt <- out.cv.glmnet$lambda.1se

# Fit the elastic net model to the data.
cat("3. Fitting elastic net model.\n")
fit.glmnet <- glmnet(X,y,family = "binomial",lambda = lambda,alpha = alpha)

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
print(with(out.cv.glmnet,
           xyplot(y ~ x,data.frame(x = log10(lambda),y = cvm),type = "l",
                  col = "blue",xlab = "log10 lambda",ylab = "error",main = "A",
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
                           type = "l",col = "blue",lty = "dashed"))),
           split = c(1,1,2,2),more = TRUE)

# Show the number of nonzero regression coefficients at different
# settings of lambda.
print(with(out.cv.glmnet,
           xyplot(y ~ x,data.frame(x = log10(lambda),y = nzero),type = "l",
                  col = "blue",main = "B",xlab = "log10 lambda",
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
               col = "blue",main = "C",xlab = "log10 lambda",
               ylab = "regression coefficient",
               scales = list(x = list(limits = c(-2.1,0.3)),
                             y = list(limits = c(-0.8,1.2))),
               panel = function(x, y, ...) {
                 panel.xyplot(x,y,...);
                 panel.abline(v = log10(lambda.opt),col = "orangered",
                              lwd = 1,lty = "dotted");
                 ltext(x = 0,y = b[nrow(b),],labels = colnames(b),pos = 4,
                       offset = 0.25,cex = 0.5)
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
fit.varbvs <- varbvs(X,NULL,y,"binomial",logodds = logodds,sa = 1,
                     verbose = FALSE)

# Compute estimates of the disease outcome using the fitted model, and
# compare against the observed values.
cat("6. Summarizing results of varbvs analysis.\n")
cat("classification results from fitted varbvs model:\n")
y.varbvs <- predict(fit.varbvs,X)
print(table(true = factor(y),pred = factor(y.varbvs)))

# Plot evolution of posterior inclusion probabilities (PIPs) at
# different settings of the prior log-odds.
trellis.device(height = 5,width = 4)
trellis.par.set(par.xlab.text = list(cex = 0.65),
                par.ylab.text = list(cex = 0.65),
                axis.text     = list(cex = 0.65))
m     <- length(logodds)
n     <- 5
vars  <- order(fit.varbvs$alpha[,m],decreasing = TRUE)[1:n]
alpha <- t(fit.varbvs$alpha[vars,])
r     <- xyplot(y ~ x,data.frame(x = logodds,y = log10(alpha[,1])),
                scales = list(x = list(limits = c(-3.6,-1.25)),
                              y = list(limits = c(-3.5,0.4))),
                type = "l",col = "blue",xlab = "prior log-odds",
                ylab = "log10 PIP",main = "D",
                panel = function(x, y, ...) {
                  panel.xyplot(x,y,...);
                  ltext(x = -1.5,y = log10(alpha[m,]),
                        labels = colnames(alpha),pos = 4,
                        offset = 0.25,cex = 0.5)
                })
for (i in 2:n)
  r <- r + as.layer(xyplot(y ~ x,
                           data.frame(x = logodds,y = log10(alpha[,i])),
                           type = "l",col = "blue"))
print(r,split = c(1,1,1,2),more = TRUE)
rm(m,n,i,vars,alpha,r)

# Show probability density of prior log-odds.
w        <- normalizelogweights(fit.varbvs$logw)
names(w) <- logodds
print(xyplot(y ~ x,data.frame(x = logodds,y = w),type = "l",col = "blue",
             xlab = "prior log-odds",ylab = "posterior prob.",main = "E",
             scales = list(x = list(limits = c(-3.6,-1.25)))) +
      as.layer(xyplot(y ~ x,data.frame(x = logodds,y = w),col = "blue",
               pch = 20,cex = 0.65)),
      split = c(1,2,1,2),more = FALSE)
rm(w)
