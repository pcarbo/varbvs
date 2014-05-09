% This is a small script to evaluate the variational approximation for the
% Bayesian variable selection mixture model in an idealized setting when the
% genetic markers (SNPs) are uncorrelated with each other.
clear

% SCRIPT PARAMETERS.
n  = 5e3;  % Number of samples.
p  = 1e3;  % Number of variables (SNPs).
na = 10;   % Number of QTLs.
r  = 0.4;  % Proportion of variance explained by "background" effects.
sd = 0.5;  % Standard deviation of additive QTL effects.

% Candidate values for the prior log-odds of inclusion (theta0), the
% proportion of variance explained by the "background" polygenic effects
% (h), and the prior standard deviation of the additive QTL effects (sa).
theta0 = (-3:0.25:-1.5)';
h      = (0.3:0.025:0.45)';
sa     = (0.2:0.1:0.8)';

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
% proportion of variance explained ("r"). That is, we adjust beta so that
%         
%   r = a/(a+1)
%
% where I've defined
%          
%   a = beta'*cov(X)*beta.
%
% Here, sb is the variance of the "background" polygenic effects.
sb   = r/(1-r)/var(X*beta,1);
beta = sqrt(sb) * beta;

% Generate the additive QTL effects.
qtl             = beta == 0;
beta(find(qtl)) = sd * randn(na,1);

% Generate the quantitative trait measurements, then take into account an
% intercept by centering the outcomes Y to have mean zero. Note that the
% columns of X are already centered.
y = X*beta + randn(n,1);
y = y - mean(y);

% COMPUTE VARIATIONAL ESTIMATES.
% Generate all candidate hyperparameter settings.
[THETA0 H SA] = ndgrid(theta0,h,sa);

% Get the variance of the "background" effects as a function of the
% proportion of variance explained.
sx = sum(var1(X));
q  = sigmoid10(-THETA0);
SB = H./(1-H)./(q*sx);

% Initialize storage for the unnormalized log-importance weights, posterior
% inclusion probabilities, and posterior means of the additive QTL effects.
ns    = numel(H);
logw  = zeros(size(H));
alpha = zeros(p,ns);
mu1   = zeros(p,ns);

% Compute the unnormalized log-importance weights. These are valid
% importance weights so long as the prior and proposal distributions cancel
% out in the expression for the importance weights.
fprintf('Computing importance weights for %d combinations ',ns);
fprintf('of hyperparameters.\n');
options.verbose = false;
fprintf('iter theta0    h   sa sigma\n');
fprintf('---- ------ ---- ---- -----\n');
for i = 1:ns
  [logw(i) alpha(:,i) mu1(:,i) mu2 s1 s2 sigma] = ...
      varbvsmix(X,y,var(y),SA(i)^2,SB(i),log(10)*THETA0(i),options);
  fprintf('%4d %0.3f %0.2f %0.2f %5.2f',i,THETA0(i),H(i),SA(i),sqrt(sigma));
  fprintf(repmat('\b',1,27))
end
fprintf('\n\n');

% Compute the normalized importance weights.
w = normalizelogweights(logw);

% Show the posterior distribution of the prior log-odds of inclusion (theta0). 
fprintf('Posterior of theta0:\n');
fprintf('theta0 prob\n')
fprintf('------ ----\n')
fprintf('%6.2f %0.2f\n',[theta0 ndsum(w,1)]')
fprintf('\n');

% Show the posterior distribution of the prior estimate of the proportion
% of variance explained by the "background" polygenic effects (h).
fprintf('Posterior of h:\n');
fprintf('  h   prob\n')
fprintf('----- ----\n')
fprintf('%0.3f %0.2f\n',[h ndsum(w,2)]')
fprintf('\n');

% Show the posterior distribution of the prior estimate of the standard
% deviation of the additive QTL effects (sa).
fprintf('Posterior of sa:\n');
fprintf(' sa  prob\n')
fprintf('---- ----\n')
fprintf('%0.2f %0.2f\n',[sa ndsum(w,3)]')
fprintf('\n');

% Calculate the posterior inclusion probabilities (PIPs), and show the
% extent to which the SNPs with high PIP correspond to the QTLs.
PIP     = alpha * w(:);
I       = find(PIP > 0.5 | qtl);
[ans J] = sort(-abs(beta(I)));
I       = I(J);
fprintf('QTL discovery:\n');
fprintf(' PIP    mu  beta\n');
fprintf('%0.2f %+0.2f %+0.2f\n',[PIP(I) mu1(I) beta(I)]');
