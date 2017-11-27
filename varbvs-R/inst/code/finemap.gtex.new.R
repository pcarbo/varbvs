# Small script to compare the fully-factorized variational
# approximation against the SuSiE ("Sum of Single Effects") approach
# for fine-mapping expression of gene FMO@ using genotype data from
# the GTEx study. The phenotype data are simulated.
library(varbvs)
library(ggplot2)
library(cowplot)
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
beta[c(201,562,952)] <- c(2,-3,4)

set.seed(1)

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
fit2 <- varbvssparse(X,y,k,sigma,sa,logodds,tol = 1e-6,maxiter = 50)

# COMPARE RESULTS
# ---------------
# Not sure it is appropriate to do so, but I thought it might be
# interesting to compare the "fit" of the two models by comparing
# their variational lower bounds (approximate marginal likelihoods).
cat("Approximate marginal log-likelihoods:\n")
last <- function (x) x[length(x)]
cat(sprintf("varbvsnorm   = %0.6f.\n",last(fit1$logw)))
cat(sprintf("varbvssparse = %0.6f.\n",last(fit2$logw)))
cat(sprintf("Approximate improvement in likelihood using SuSie model: %0.2f\n",
            exp(last(fit2$logw) - last(fit1$logw))))

# Show convergence of the co-ordinate ascent updates for the two methods.
p1 <- ggplot(data.frame(x = 1:length(fit1$logw),
                        y = log10(max(fit1$logw) - fit1$logw)),
             aes(x = x,y = y)) +
    geom_line(color = "dodgerblue",size = 1) +
    labs(x = "iteration",y = "log10-distance from max",
         title = "fully-factorized")
p2 <- ggplot(data.frame(x = 1:length(fit2$logw),
                        y = log10(max(fit2$logw) - fit2$logw)),
             aes(x = x,y = y)) +
    geom_line(color = "darkorange",size = 1) +
    labs(x = "iteration",y = "",title = "sum of single effects")
plot_grid(p1,p2)

# Plot the posterior inclusion probabilities (PIPs) under the
# fully-factorized variational approximation.
clrs    <-   c("#0D0887FF","#5402A3FF","#8B0AA5FF","#B93289FF",
               "#DB5C68FF","#F48849FF","#FEBC2AFF","#F0F921FF")
markers <- which(beta != 0)
r2 <- apply((cor(X)^2)[,markers],1,max)
p3 <- ggplot(data.frame(x = 1:p,y = fit1$alpha,r2 = r2),
             aes(x = x,y = y,color = r2)) +
    geom_point(shape = 20,size = 2) +
    scale_x_continuous(breaks = NULL) +
    scale_color_gradientn(colors = clrs) +
    labs(x = "",y = "PIP",title = "fully-factorized")

# Plot the posterior inclusion probabilities (PIPs) under the SuSie
# model. We report the PIPs as the probability that at least one of
# the coefficients is nonzero.
pip <- 1 - apply(1 - fit2$alpha,1,prod)
p4 <- ggplot(data.frame(x = 1:p,y = pip,r2 = r2),
             aes(x = x,y = y,color = r2)) +
    geom_point(shape = 20,size = 2) +
    scale_x_continuous(breaks = NULL) +
    scale_color_gradientn(colors = clrs) +
    labs(x = "",y = "PIP",title = "sum of single effects")

print(plot_grid(p3,p4,nrow = 2))

sessionInfo()

