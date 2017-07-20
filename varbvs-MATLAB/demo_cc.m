% This script illustrates 'varbvs' for genome-wide mapping of a binary
% (e.g., case-control) trait in a simulated data set in which all the
% genetic markers are uncorrelated with each other (i.e., they are
% "unlinked").
%
% Part of the varbvs package, https://github.com/pcarbo/varbvs
%
% Copyright (C) 2012-2017, Peter Carbonetto
%
% This program is free software: you can redistribute it under the
% terms of the GNU General Public License; either version 3 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANY; without even the implied warranty of
% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
%
clear

% SCRIPT PARAMETERS
% -----------------
n  = 1500;  % Number of samples (subjects).
p  = 2000;  % Number of variables (genetic markers).
m  = 2;     % Number of covariates (m >= 0).
na = 20;    % Number of markers that affect the binary outcome.
sa = 0.15;  % Variance of log-odds ratios.
p1 = 0.25;  % Target proportion of subjects that are cases (y = 1).

% Candidate values for the prior log-odds of inclusion.
logodds = (-3:0.1:-1)';

% Set the random number generator seed.
rng(1);

% GENERATE DATA SET
% -----------------
% Generate the minor allele frequencies so that they are uniform over range
% [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
% (i.e., unlinked), according to the specified minor allele frequencies.
fprintf('1. GENERATING DATA SET.\n')
maf = 0.05 + 0.45 * rand(1,p);
X   = (rand(n,p) < repmat(maf,n,1)) + ...
      (rand(n,p) < repmat(maf,n,1));
X   = single(X);

% Generate additive effects for the markers so that exactly na of them
% have a nonzero effect on the trait.
i       = randperm(p);
i       = i(1:na);
beta    = zeros(p,1);
beta(i) = sqrt(sa)*randn(na,1);

% Generate random labels for the markers.
labels = num2cell(randi(max(1e6,p),p,1));
labels = cellfun(@num2str,labels,'UniformOutput',false);
labels = strcat('rs',labels);

% Generate the covariate data (Z), and the linear effects of the
% covariates (u).
if m > 0
  Z = randn(n,m);
  u = randn(m,1);
else
  Z = [];
end

% For each sample, calculate the probability of being a case (y = 1).
mu = logit(p1);
t  = mu + X*beta;
if m > 0
  t = t + Z*u;
end
t = double(t);

% Simulate the binary trait (case-control status) as a coin toss with
% success rates given by the logistic regression.
y = rand(n,1) < sigmoid(t);

% FIT VARIATIONAL APPROXIMATION TO POSTERIOR
% ------------------------------------------
% Fit the fully-factorized variational approximation to the posterior
% distribution of the coefficients for a logistic regression model of a
% binary outcome (case-control status), with spike and slab priors on the
% coefficients.
fprintf('2. FITTING MODEL TO DATA.\n')
fit = varbvs(X,Z,y,labels,'binomial',struct('logodds',logodds));
fprintf('\n');

% SUMMARIZE POSTERIOR DISTRIBUTION
% --------------------------------
fprintf('3. SUMMARIZING RESULTS.\n')
varbvsprint(fit);
fprintf('\n');

% EVALUATE MODEL PREDICTIONS
% --------------------------
fprintf('4. EVALUATING FITTED MODEL.\n');
ypred = varbvspredict(fit,X,Z);
fprintf('     predicted\n');
fprintf('true    0    1\n');
fprintf('  0  %4d %4d\n',sum(y == 0 & ypred == 0),sum(y == 0 & ypred == 1));
fprintf('  1  %4d %4d\n',sum(y == 1 & ypred == 0),sum(y == 1 & ypred == 1));
