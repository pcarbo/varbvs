# In this small example, we explore the posterior distribution of the
# coefficients in a linear regression, with spike and slab priors. In this
# idealized case, the variables are independent.
source("varbvs.R")

# SCRIPT PARAMETERS.
p  <- 1e3  # The number of variables (SNPs).
n  <- 500  # The number of samples.
na <- 20   # Number of variables that effect the outcome ("causal" SNPs).

# The level of residual noise is set so that the proportion of variance
# explained should be a little less than 30% on average (according to my
# rough calculations).
se <- 9

# These two parameters specify the Beta prior on the proportion of variables
# (SNPs) that are included in the linear model of Y.
a <- 0.02
b <- 1

# This parameter specifies the prior on the variance of the regression
# coefficients (sa). For more information on this prior, see MULTISNPSIM.
ca <- 0.02

# Candidate values of the variance of the residual (sigma), the prior
# variance of the regression coefficients (sa), and logarithm of the prior
# inclusion probability (log10q).
sigma  <- seq(8,13)
sa     <- seq(0.025,0.4,0.025)
log10q <- seq(-2.5,-1,0.25)

# Set the random number generator seed.
set.seed(1);

# CREATE THE DATA.
# Note that X and y are centered.
print(noquote("Creating data."))
snps <- create.snps(p,na)
data <- create.data(snps,se,n)
