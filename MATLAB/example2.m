% In this small example, I simulate a data set from an idealized genetic
% association study in which the genetic markers (single nucleotide
% polymorphisms, or SNPs) are uncorrelated (i.e. no linkage disequilibrium),
% and I assess the performance of the variational approximation for mapping
% associations between the genetic markers and a quantitative trait in
% this setting.
%
% Note that I have tested this script in MATLAB 8.1, and it may not work
% correctly in earlier versions of MATLAB.
clear

% SCRIPT PARAMETERS.
p  = 1000;  % Number of variables (SNPs).
n  = 2e3;   % Number of samples (subjects).
na = 10;    % Number of SNPs that affect the outcome.
sb = 0.2;   % Standard deviation of nonzero coefficients.
se = 9;     % Variance of residual.

% Candidate values for the prior proportion of variance explained (h), and
% prior log-odds for inclusion (theta0).
h      = (0.1:0.1:0.5)';
theta0 = (-2.5:0.25:-1.5)';

% Set the random number generator seed.
seed = 1;
rng(seed);

% CREATE THE DATA.
% Note that X and y are centered.
fprintf('Creating data.\n');
[maf beta] = createsnps(p,na);
beta       = sb * beta;
[X y]      = createdata(maf,beta,se,n);
