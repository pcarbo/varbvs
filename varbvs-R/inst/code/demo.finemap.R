# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(cowplot)
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
# This vector specifies the additive effects of the SNPs.
beta            <- rep(0,531)
beta[c(88,506)] <- c(0.25,-0.35)

# LOAD GENOTYPE DATA
# ------------------
# Load the genotype data for the 531 SNPs at the bone-mineral density
# locus on chromosome 11 that was identified in the Nature Genetics
# manuscript.
cat("Loading genotype data.\n")
load("../datafiles/cfw.chr11bmd.RData")
n <- nrow(X)
p <- ncol(X)

# Center the columns of X so that each column has a mean of zero.
rep.row <- function (x, n)
  matrix(x,n,length(x),byrow = TRUE)
X <- X - rep.row(colMeans(X),n)

# Set the random number generator seed.
set.seed(1)

# SIMULATE PHENOTYPE DATA
# -----------------------
# Generate the quantitative trait measurements.
cat("Simulating phenotype data.\n")
y <- c(X %*% beta + rnorm(n))

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of a
# continuous outcome (quantitiative trait), with spike and slab priors on
# the coefficients. 
cat("Fitting model to data.\n")
fit <- varbvs(X,NULL,y,sa = 1,logodds = log10(1/p),verbose = TRUE)

# COMPUTE "PROXY" POSTERIOR PROBABILITIES
# ---------------------------------------
sigmoid <- function (x)
  1/(1 + exp(-x))

# TO DO: Explain what this code chunk does.
top.markers <- summary(fit)$top.vars$index
i <- top.markers[1]
j <- top.markers[2]

# Compute a few useful quantities. 
xy <- c(y %*% X)
d  <- diag(t(X) %*% X)

# TO DO: Explain what this code chunk does.
pip1    <- c(fit$alpha)
pip1[i] <- 0
Xr      <- c(X %*% (pip1 * c(fit$mu)))
BF      <- rep(0,p)
    
# Repeat for each candidate marker.
for (k in 1:p) {
    
  # Update the variational estimate of the posterior mean.
  b <- with(fit,s[k]/sigma * (xy[k] - sum(X[,k]*Xr)))
        
  # Compute the variational estimate of the posterior inclusion
  # probability.
  BF[k] <- with(fit,s[k]/(sigma*sa) * exp(b^2/(2*s[k])))
}

ppp <- BF / sum(BF)

# PLOT RESULTS
# ------------
# TO DO: Explain what this code chunk does.
cat("Plotting results.\n")
R <- cor(X)

# TO DO: Explain what this plot shows.
clrs <- c("midnightblue","darkviolet","darkorchid","maroon","tomato","orange")
regionplot1 <- ggplot(data.frame(pos = 1:p,pip = fit$pip),
                     aes(x = pos,y = pip,color = abs(R[,i]))) +
  geom_point(shape = 20,size = 2.5) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(trans = "log10",breaks = 10^(-4:0)) +
  scale_color_gradientn(colors = clrs) +
  labs(x = "position on chromosome 11, 94-98 Mb",y = "varbvs PIP",
       title = paste("correlations with",map[i,"id"])) +
  theme_cowplot(font_size = 11) +
  theme(axis.line = element_blank())

# TO DO: Explain what this plot shows.
regionplot2 <- ggplot(data.frame(pos = 1:p,pip = fit$pip),
                     aes(x = pos,y = pip,color = abs(R[,j]))) +
  geom_point(shape = 20,size = 2.5) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(trans = "log10",breaks = 10^(-4:0)) +
  scale_color_gradientn(colors = clrs) +
  labs(x = "position on chromosome 11, 94-98 Mb",y = "varbvs PIP",
       title = paste("correlations with",map[j,"id"])) +
  theme_cowplot(font_size = 11) +
  theme(axis.line = element_blank())

# Arrange the plots.
print(plot_grid(regionplot1,regionplot2,nrow = 2))

