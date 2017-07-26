% Functions "varbvs" and "varbvsmix" should produce the same estimates when
% there are exactly two mixture components (a "spike" and a "slab"). This
% script verifies this in a simulated data set.
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

% SCRIPT PARAMETERS
% -----------------
n  = 1000;  % Number of samples.
m  = 2;     % Number of covariates (m >= 0).
p  = 2000;  % Number of variables (genetic markers).
se = 4;     % Variance of residual.

% The standard deviations and mixture weights used to simulate the additive
% effects on the quantitative trait. Note that the first mixture component
% must have a standard deviation of exactly zero.
sd = [    0  0.5 ]';
q  = [ 0.95 0.05 ]';

% Set the random number generator seed.
rng(1);

% GENERATE DATA SET
% -----------------
% Generate the minor allele frequencies so that they are uniform over range
% [0.05,0.5]. Then simulate genotypes assuming all markers are uncorrelated
% (i.e., unlinked), according to the specified minor allele frequencies.
fprintf('1. GENERATING DATA SET.\n');
fprintf('Data simulation settings:\n');
fprintf('  - Num. data samples       %d\n',n);
fprintf('  - Num. covariates         %d\n',m);
fprintf('  - Num. variables (SNPs)   %d\n',p);
fprintf('  - Num. mixture components %d\n',length(q));
maf = 0.05 + 0.45 * rand(1,p);
X   = (rand(n,p) < repmat(maf,n,1)) + ...
      (rand(n,p) < repmat(maf,n,1));
X   = single(X);

% Generate additive effects according to the specified standard deviations
% (sd) and mixture weights (q).
k    = randtable(q,p);
beta = sd(k) .* randn(p,1);

% Generate random labels for the markers.
labels = num2cell(randi(max(1e6,p),p,1));
labels = cellfun(@num2str,labels,'UniformOutput',false);
labels = strcat('rs',labels);

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
% Fit the varbvsmix model.
fprintf('2. FITTING VARBVSMIX MODEL TO DATA.\n');
alpha   = rand(p,1);
alpha   = [ alpha 1-alpha ];
mu      = [ zeros(p,1) randn(p,1) ];
options = struct('alpha',alpha,'mu',mu);
fit     = varbvsmix(X,Z,y,sd.^2,labels,options);

% Fit the varbvs model.
fprintf('3. FITTING VARBVS MODEL TO DATA.\n');
q1      = fit.q(2);
options = struct('sa',sd(2)^2,'logodds',log10(q1/(1-q1)),...
                 'alpha',alpha(:,2),'mu',mu(:,2));
fit2    = varbvs(X,Z,y,labels,[],options);

% Check that the parameter estimates are the same.
relerr = @(x, y) abs(x - y)/min(abs(x),abs(y));
fprintf('Max. difference in PIPs:             %0.2e\n',...
        max(abs(fit.alpha(:,2) - fit2.alpha)));
fprintf('Max. difference in posterior means:  %0.2e\n',...
        max(abs(fit.mu(:,2) - fit2.mu)));
fprintf('Relative difference in lower bounds: %0.2e\n',...
        relerr(fit.logw(end),fit2.logw));