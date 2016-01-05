% This is a small example to check that estimation of the residual
% variance parameter (SIGMA) works correctly.
clear

% SCRIPT PARAMETERS.
p  = 1e3;  % The number of variables (SNPs).
n  = 2e3;  % The number of samples.
na = 20;   % Number of variables that affect the outcome.
sa = 0.1;  % Prior variance of the regression coefficients.

% The level of residual noise is set so that the proportion of variance
% explained should be a little less than 30% on average (according to my
% rough calculations).
se = 9;

% Set the random number generator seed.
seed = 1;
rand('state',seed);
randn('state',seed);

% CREATE THE DATA.
% Note that X and y are centered.
fprintf('Creating data.\n');
[maf beta] = createsnps(p,na);
[X y]      = createdata(maf,beta,se,n);

% COMPUTE VARIATIONAL ESTIMATES.
fprintf('Computing variational estimates.\n');
options.update_sigma = true;
[lnZ alpha mu s sigma] = varbvs(X,y,se,sa,log10(na/p),options);
