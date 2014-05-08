% This is a small script to evaluate the variational approximation for the
% Bayesian variable selection mixture model in an idealized setting when the
% genetic markers (SNPs) are uncorrelated with each other.
clear

% SCRIPT PARAMETERS.
n  = 1e3;  % Number of samples.
p  = 2e3;  % Number of variables (SNPs).
na = 10;   % Number of QTLs.
h  = 0.4;  % Proportion of variance explained by "background" effects.
sa = 0.2;  % Variance of additive QTL effects.

% Set the random number generator seed.
seed = 1;
rand('state',seed);
randn('state',seed);

% GENERATE DATA SET.
% Generate the minor allele frequencies (MAFs) so that they are uniform over
% range [0.05,0.5]. Simulate genotype data X assuming all the SNPs are
% uncorrelated (i.e. no linkage disequilibrium between the SNPs), according
% to the specified MAFs, then center the columns of X.
fprintf('Creating data set.\n');
maf = 0.05 + 0.45 * rand(1,p);
X   = (rand(n,p) < repmat(maf,n,1)) + ...
      (rand(n,p) < repmat(maf,n,1));
X   = X - repmat(mean(X),n,1);

% Generate "background" polygenic effects for most (but not all) of the
% SNPs.
beta    = randn(p,1);
I       = randperm(p);
I       = I(1:na);
beta(I) = 0;

% Adjust the "background" polygenic effects so that we control for the
% proportion of variance explained ("h"). That is, we adjust beta so that
%         
%   h = a/(a+1)
%
% where I've defined
%          
%   a = beta'*cov(X)*beta.
%
% Here, sb is the variance of the "background" polygenic effects.
sb   = h/(1-h)/var(X*beta,1);
beta = sqrt(sb) * beta;

% Generate the additive QTL effects.
qtl     = beta == 0;
I       = find(qtl);
beta(I) = sqrt(sa) * randn(na,1);

% Generate the quantitative trait measurements, then take into account an
% intercept by centering the outcomes Y to have mean zero. Note that the
% columns of X are already centered.
y = X*beta + randn(n,1);
y = y - mean(y);

% COMPUTE VARIATIONAL ESTIMATES.
fprintf('Computing variational estimates.\n');
[lnZ alpha mu1 mu2 s1 s2] = varbvsmix(X,y,var(y),sa,sb,log(na/p));
