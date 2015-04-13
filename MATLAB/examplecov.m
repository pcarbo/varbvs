% EXPLAIN HERE WHAT THIS SCRIPT DOES.
%
% This script has been tested in MATLAB R2014b (8.4).
%
clear

% SCRIPT PARAMETERS
% -----------------
n  = 1e3;  % Number of samples.
p  = 2e3;  % Number of variables (markers).
m  = 3;    % Number of covariates (variables included in model with Pr. 1).
na = 20;   % Number of QTLs.
r  = 0.5;  % Proportion of variance in residuals explained by QTLs.

% Candidate values for the prior proportion of variance explained (h) and
% prior log-odds of inclusion (theta0).
theta0 = (-3:0.5:-1)';
h      = (0.1:0.1:0.9)';

% Set the random number generator seed.
rng(1);

% GENERATE DATA SET
% -----------------
% Generate the minor allele frequencies so that they are uniform over range
% [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
% (i.e. no linkage disequilibrium), according to the specified minor allele
% frequencies.
fprintf('Creating data set.\n');
maf = 0.05 + 0.45 * rand(1,p);
X   = (rand(n,p) < repmat(maf,n,1)) + ...
      (rand(n,p) < repmat(maf,n,1));

% Generate the covariates, and the linear effects of the covariates.
Z = randn(n,m);
u = randn(m,1);

% Generate additive effects for the markers, such that 'na' of them have a
% nonzero effect on the trait.
I       = randperm(p);
I       = I(1:na);
beta    = zeros(p,1);
beta(I) = randn(na,1);

% Adjust the QTL effects so that we control for the proportion of variance
% explained (r). That is, we adjust beta so that
%         
%   r = a/(a+1)
%
% where I've defined
%          
%   a = beta'*cov(X)*beta.
%
% Here, sb is the variance of the (nonzero) QTL effects.
sb   = r/(1-r)/var(X*beta,1);
beta = sqrt(sb) * beta;

% Generate the quantitative trait measurements.
y = Z*u + X*beta + randn(n,1);

% Adjust the genotypes and phenotypes so that the linear effects of
% the covariates are removed. This is equivalent to integrating out
% the regression coefficients corresponding to the covariates with
% respect to an improper, uniform prior; see Chipman, George and
% McCulloch, "The Practical Implementation of Bayesian Model
% Selection," 2001.
%
% Note that this should give the same result as centering the
% columns of X and subtracting the mean from y when we have only
% one covariate, the intercept.
y = y - Z*((Z'*Z)\(Z'*y));
X = X - Z*((Z'*Z)\(Z'*X));

% COMPUTE VARIATIONAL APPROXIMATION
% ---------------------------------
fprintf('Computing variational estimates.\n');
[H THETA0] = ndgrid(h,theta0);
[logw sigma alpha mu s] = multisnphyper(X,y,H,THETA0,false);

% Compute the normalized importance weights.
w = normalizelogweights(logw);

% SUMMARIZE RESULTS OF MULTI-MARKER ANALYSIS
% ------------------------------------------
% Show the posterior mean of each of the hyperparameters.
fprintf('Posterior mean of hyperparameters:\n');
fprintf('param   mean\n');
fprintf('sigma  %5.2f\n',dot(sigma(:),w(:)));
fprintf('theta0 %+0.2f\n',dot(THETA0(:),w(:)));
fprintf('h      %5.2f\n',dot(H(:),w(:)));
fprintf('\n');

% Calculate the posterior inclusion probabilities (PIPs), and show the PIPs
% for the variables that are most likely to be included in the model,
% with PIP > 0.1.
PIP     = alpha * w(:);
[ans I] = sort(-PIP);
I       = I(1:sum(PIP > 0.1));
fprintf('Selected variables:\n');
fprintf(' PIP    mu  beta\n');
fprintf('%0.2f %+0.2f %+0.2f\n',[PIP(I) mu(I) beta(I)]');
