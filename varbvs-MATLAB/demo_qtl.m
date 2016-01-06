% This script illustrates 'varbvs' genome-wide mapping of a quantitative
% trait in simulated data when all the genetic markers are "unlinked" (not
% correlated with each other).
clear

% SCRIPT PARAMETERS
% -----------------
n  = 500;  % Number of samples.
p  = 2e3;  % Number of variables (genetic markers).
m  = 3;    % Number of covariates (m >= 0).
na = 20;   % Number of quantitative trait loci (QTLs).
se = 4;    % Variance of residual.
r  = 0.5;  % Proportion of variance in trait explained by QTLs.

% Candidate values for the prior log-odds of inclusion.
theta0 = (-3:0.1:-1.5)';

% Set the random number generator seed.
rng(1);

% GENERATE DATA SET
% -----------------
% Generate the minor allele frequencies so that they are uniform over range
% [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
% (i.e. no linkage disequilibrium), according to the specified minor allele
% frequencies.
fprintf('Generating data set.\n');
maf = 0.05 + 0.45 * rand(1,p);
X   = (rand(n,p) < repmat(maf,n,1)) + ...
      (rand(n,p) < repmat(maf,n,1));

% Generate additive effects for the markers so that exactly na of them have
% a nonzero effect on the trait.
I       = randperm(p);
I       = I(1:na);
beta    = zeros(p,1);
beta(I) = randn(na,1);

% Adjust the QTL effects so that we control for the proportion of variance
% explained (r). That is, we adjust beta so that r = a/(a+1), where I've
% defined a = beta'*cov(X)*beta. Here, sb is the variance of the (nonzero)
% QTL effects.
sb   = r/(1-r)/var(X*beta,1);
beta = sqrt(sb*se) * beta;

% Generate the covariate data (Z), and the linear effects of the
% covariates (u).
if m > 0
  Z = randn(n,m);
  u = randn(m,1);
else
  Z = [];
end
  
% Generate the quantitative trait measurements.
y = X*beta + sqrt(se)*randn(n,1);
if m > 0
  y = y + Z*u;
end
