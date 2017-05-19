# The varbvs and varbvsmix functions should produce the same estimates
# when there are exactly 2 mixture components (a "spike" and a
# "slab"). This script verifies this in a small simulated data set.

# SCRIPT PARAMETERS
# -----------------
n  = 1000;  # Number of samples.
m  = 2;     # Number of covariates (m >= 0).
p  = 2000;  # Number of variables (genetic markers).
se = 4;     # Variance of residual.

% The standard deviations and mixture weights used to simulate the additive
% effects on the quantitative trait. Note that the first mixture component
% must have a standard deviation of exactly zero.
sd = [    0  0.5 ]';
q  = [ 0.95 0.05 ]';

% Set the random number generator seed.
rng(1);
