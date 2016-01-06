% This script illustrates 'varbvs' for genome-wide mapping of a quantitative
% trait in a simulated data set in which all the genetic markers are
% uncorrelated with each other (i.e., they are "unlinked").
clear

% SCRIPT PARAMETERS
% -----------------
n  = 500;  % Number of samples.
p  = 2e3;  % Number of variables (genetic markers).
m  = 3;    % Number of covariates (m >= 0).
na = 20;   % Number of quantitative trait loci (QTLs).
se = 4;    % Variance of residual.
r  = 0.5;  % Proportion of variance in trait explained by QTLs.

% Include an intercept?
intercept = true;

% Candidate values for the prior log-odds of inclusion.
logodds = (-3:0.1:-1)';

% Set the random number generator seed.
rng(1);

% GENERATE DATA SET
% -----------------
% Generate the minor allele frequencies so that they are uniform over range
% [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
% (i.e., unlinked), according to the specified minor allele frequencies.
maf = 0.05 + 0.45 * rand(1,p);
X   = (rand(n,p) < repmat(maf,n,1)) + ...
      (rand(n,p) < repmat(maf,n,1));

% Generate additive effects for the markers so that exactly na of them have
% a nonzero effect on the trait.
i       = randperm(p);
i       = i(1:na);
beta    = zeros(p,1);
beta(i) = randn(na,1);

% Adjust the QTL effects so that we control for the proportion of variance
% explained (r). That is, we adjust beta so that r = a/(a+1), where I've
% defined a = beta'*cov(X)*beta. Here, sb is the variance of the (nonzero)
% QTL effects.
sb   = r/(1-r)/var(X*beta,1);
beta = sqrt(sb*se) * beta;

% Generate the intercept.
if intercept
  mu = randn;
end

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
if intercept
  y = y + mu;
end
if m > 0
  y = y + Z*u;
end

% FIT VARIATIONAL APPROXIMATION TO POSTERIOR
% ------------------------------------------
% Fit the fully-factorized variational approximation to the posterior
% distribution of the coefficients for a linear regression model of a
% continuous outcome (quantitiative trait), with spike and slab priors on
% the coefficients.
fit = varbvs(X,Z,y,[],struct('intercept',intercept,'logodds',logodds));
