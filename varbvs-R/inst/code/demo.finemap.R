# TO DO: Explain here what this script does, and how to use it.
library(ggplot2)
library(cowplot)
library(varbvs)

# SCRIPT PARAMETERS
# -----------------
# This vector specifies the additive effects of the SNPs.
beta             <- rep(0,531)
beta[c(161,506)] <- c(0.4,-0.4)

# Set the random number generator seed.
set.seed(1)

# LOAD GENOTYPE DATA
# ------------------
# Load the genotype data for the 531 SNPs at the bone-mineral density
# locus on chromosome 11 that was identified in the Nature Genetics
# manuscript.
cat("Loading genotype data.\n")
load("../datafiles/cfw.chr11bmd.RData")

# SIMULATE PHENOTYPE DATA
# -----------------------
# Generate the quantitative trait measurements.
cat("Simulating phenotype data.\n")
n <- nrow(X)
p <- ncol(X)
y <- c(X %*% beta + rnorm(n))

# FIT VARIATIONAL APPROXIMATION TO POSTERIOR
# ------------------------------------------
# Fit the fully-factorized variational approximation to the posterior
# distribution of the coefficients for a linear regression model of a
# continuous outcome (quantitiative trait), with spike and slab priors on
# the coefficients. 
cat("Fitting model to data.\n")
fit <- varbvs(X,NULL,y,sa = 1,logodds = log10(1/p),verbose = FALSE)

# PLOT RESULTS
# ------------
# TO DO: Explain what this plot shows.
R    <- cor(X)
i    <- order(fit$pip,decreasing = TRUE)
r2   <- R[,i[1]]^2
clrs <- c("midnightblue","darkviolet","darkorchid","maroon","tomato","orange")
regionplot1 <- ggplot(data.frame(pos = 1:p,pip = fit$pip),
                     aes(x = pos,y = pip,color = r2)) +
  geom_point(shape = 20,size = 2.5) +
  scale_x_continuous(breaks = NULL) +
  scale_y_continuous(trans = "log10",breaks = 10^(-4:0)) +
  scale_color_gradientn(colors = clrs) +
  labs(x = "position on chromosome 11, 94-98 Mb",y = "PIP") +
  theme_cowplot(font_size = 12) +
  theme(axis.line = element_blank())
regionplot1
