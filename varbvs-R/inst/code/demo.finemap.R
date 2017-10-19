# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(cowplot)
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
# This vector specifies the additive effects of the SNPs.
beta            <- rep(0,531)
beta[c(88,506)] <- c(0.2,-0.35)

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
fit <- varbvs(X,NULL,y,sa = 1,verbose = FALSE)

# COMPUTE PROXY PROBABILITIES
# ---------------------------
# TO DO: Explain what this code chunk does.
cat("Computing proxy probabilities for the top 2 SNPs.\n")
source("finemap.R")
top.markers <- summary(fit)$top.vars$index
i     <- top.markers[1]
j     <- top.markers[2]
vars1 <- which(R[,i] >= 0.1)
vars2 <- which(R[,j] >= 0.1)
BF1   <- varbvsproxybf(X,NULL,y,fit,i,vars1) 
BF2   <- varbvsproxybf(X,NULL,y,fit,j,vars2)

# TO DO: Explain what this code chunk does.
bf1 <- BF1[,1]
bf2 <- BF2[,1]
pp1 <- bf1/sum(bf1)
pp2 <- bf2/sum(bf2)

# TO DO: Explain what this code chunk does.
markers1 <- order(bf1,decreasing = TRUE)
n1       <- min(which(cumsum(pp1[markers1]) > 0.95))
markers1 <- markers1[1:n1]

# TO DO: Explain what this code chunk does.
markers2 <- order(bf2,decreasing = TRUE)
n2       <- min(which(cumsum(pp2[markers2]) > 0.95))
markers2 <- markers2[1:n2]

# PLOT RESULTS
# ------------
# TO DO: Explain what this code chunk does.
cat("Plotting results.\n")
R <- cor(X)

# TO DO: Explain what this plot shows.
clrs <- c("midnightblue","darkviolet","darkorchid","maroon",
          "tomato","orange","gold","yellow")
dat1 <- data.frame(pos = 1:p,pip = fit$pip)
regionplot1 <- ggplot(mapping = aes(x = pos,y = pip,color = abs(R[,i]))) +
  geom_point(data = dat1,size = 2,shape = 19) +
  geom_point(data = dat1[vars1[markers1],],color = "black",shape = 4,
             size = 0.5) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(trans = "log10",breaks = 10^(-4:0)) +
  scale_color_gradientn(colors = clrs) +
  labs(x = "position on chromosome 11, 94-98 Mb",y = "varbvs PIP",
       title = paste("correlations with",map[i,"id"])) +
  theme_cowplot(font_size = 11) +
  theme(axis.line = element_blank())

# TO DO: Explain what this plot shows.
dat2 <- data.frame(pos = 1:p,pip = fit$pip)
regionplot2 <- ggplot(mapping = aes(x = pos,y = pip,color = abs(R[,j]))) +
  geom_point(data = dat2,size = 2,shape = 19) +
  geom_point(data = dat2[vars2[markers2],],color = "black",shape = 4,
             size = 0.5) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(trans = "log10",breaks = 10^(-4:0)) +
  scale_color_gradientn(colors = clrs) +
  labs(x = "position on chromosome 11, 94-98 Mb",y = "varbvs PIP",
       title = paste("correlations with",map[i,"id"])) +
  theme_cowplot(font_size = 11) +
  theme(axis.line = element_blank())

# Arrange the plots.
print(plot_grid(regionplot1,regionplot2,nrow = 2))

