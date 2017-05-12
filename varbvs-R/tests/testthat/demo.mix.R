# This script illustrates the "varbvsmix" function on a simulated data
# set, in which all candidate variables (predictors) are uncorrelated.
library(varbvs)  

# SCRIPT PARAMETERS
# -----------------
n  <- 1000  # Number of samples.
m  <- 2     # Number of covariates (m >= 0).
p  <- 2000  # Number of variables (genetic markers).
se <- 4     # Variance of residual.

# The standard deviations and mixture weights used to simulate the
# additive effects on the quantitative trait. Note that the first
# mixture component must have a standard deviation of exactly zero.
sd = [    0  0.1  0.2  0.5 ]';
q  = [ 0.95 0.03 0.01 0.01 ]';

# Set the random number generator seed.
set.seed(1)

