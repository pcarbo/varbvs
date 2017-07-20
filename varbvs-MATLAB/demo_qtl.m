% This script illustrates usage of the function varbvs for genome-wide
% mapping of a quantitative trait. The data set is simulated assuming that
% all the genetic markers are uncorrelated with each other (i.e., they are
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
n  = 800;   % Number of samples.
p  = 2000;  % Number of variables (genetic markers).
m  = 3;     % Number of covariates (m >= 0).
na = 20;    % Number of quantitative trait loci (QTLs).
se = 4;     % Variance of residual.
r  = 0.5;   % Proportion of variance in trait explained by QTLs.

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

% Generate additive effects for the markers so that exactly na of them have
% a nonzero effect on the trait.
i       = randperm(p);
i       = i(1:na);
beta    = zeros(p,1);
beta(i) = randn(na,1);

% Generate random labels for the markers.
labels = num2cell(randi(max(1e6,p),p,1));
labels = cellfun(@num2str,labels,'UniformOutput',false);
labels = strcat('rs',labels);

% Adjust the QTL effects so that we control for the proportion of variance
% explained (r). That is, we adjust beta so that r = a/(a+1), where I've
% defined a = beta'*cov(X)*beta. Here, sb is the variance of the (nonzero)
% QTL effects.
sb   = double(r/(1-r)/var(X*beta,1));
beta = sqrt(sb*se) * beta;

% Generate the intercept.
mu = randn;

% Generate the covariate data (Z), and the linear effects of the
% covariates (u).
if m > 0
  Z = randn(n,m);
  u = randn(m,1);
else
  Z = [];
end
  
% Generate the quantitative trait measurements.
y = mu + X*beta + sqrt(se)*randn(n,1);
if m > 0
  y = y + Z*u;
end
y = double(y);

% FIT VARIATIONAL APPROXIMATION TO POSTERIOR
% ------------------------------------------
% Fit the fully-factorized variational approximation to the posterior
% distribution of the coefficients for a linear regression model of a
% continuous outcome (quantitiative trait), with spike and slab priors on
% the coefficients.
fprintf('2. FITTING MODEL TO DATA.\n')
fit = varbvs(X,Z,y,labels,[],struct('logodds',logodds));
fprintf('\n');

% SUMMARIZE POSTERIOR DISTRIBUTION
% --------------------------------
fprintf('3. SUMMARIZING RESULTS.\n')
gvarbvsprint(fit);
fprintf('\n');

% EVALUATE MODEL PREDICTIONS
% --------------------------
fprintf('4. EVALUATING FITTED MODEL.\n');
ypred = varbvspredict(fit,X,Z);
r     = corrcoef(y,ypred);
r     = r(2)^2;
fprintf('r^2 between predicted Y and observed Y is %0.3f.\n',r);
plot(y,ypred,'o','MarkerFaceColor',[0.10 0.10 0.44],'MarkerEdgeColor',...
     'none','MarkerSize',6)
xlabel('observed Y');
ylabel('estimated Y')
set(gca,'FontSize',12);
