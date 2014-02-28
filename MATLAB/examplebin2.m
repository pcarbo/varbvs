% This script is similar to examplebin.m, except that I now also include
% additional covariates in the logistic regression.
clear

% SCRIPT PARAMETERS.
p  = 1000;  % Number of variables (SNPs).
n  = 2e4;   % Number of samples (subjects).
na = 10;    % Number of SNPs that affect the outcome.
sb = 0.2;   % Standard deviation of log-odds ratios.
p1 = 0.4;   % Target proportion of subjects that are cases (y = 1).

% Regression coefficients corresponding to covariates (other than the
% intercept).
u = [-1 1]';

% Candidate values for the prior proportion of variance explained (h), and
% prior log-odds for inclusion (theta0).
h      = (0.1:0.1:0.5)';
theta0 = (-3:0.5:-1)';

% Set the random number generator seed.
seed = 1;
rng(seed);

% GENERATE DATA SET.
% To help ensure that the proportion of subjects with the disease ("cases")
% is equal to the target proportion (p1), I adjust the coefficients in the
% logistic regression so that every SNP that has an effect on disease status
% is matched with a SNP that has the same effect, but in the opposite
% direction. I also adjust the minor allele frequencies so that these
% "matched" SNPs have the same minor allele frequency.
fprintf('Creating data set.\n');
[maf beta]  = createsnps(p,na);
beta        = sb * beta;
I           = find(beta ~= 0);
I           = I(randperm(na));
t           = na/2;
snps1       = I(1:t);
snps2       = I(t+1:na);
beta(snps2) = -beta(snps1);
maf(snps2)  = maf(snps1);

% Generate the samples.
nz    = length(u);
[X y] = createbindata([maf repmat(0.5,1,nz)],[beta; u],logit(p1),n);
Z     = X(:,p+(1:nz));
X     = X(:,1:p);

% COMPUTE VARIATIONAL ESTIMATES.
fprintf('Computing variational estimates.\n');
[H THETA0] = ndgrid(h,theta0);
[logw alpha mu s eta] = multisnpbinzhyper(X,Z,y,H,THETA0);

% Compute the normalize importance weights.
w = normalizelogweights(logw);

% Calculate the posterior mean of the prior variance of the log-odds
% ratios (sa).
sx = sum(var1(X));
sa = pve2sa(sx,H,THETA0);
fprintf('Posterior mean of sa: %0.2f\n\n',dot(sa(:),w(:)));

% Show the posterior distribution of the genome-wide log-odds (theta0). 
fprintf('Posterior of theta0:\n');
fprintf('theta0 prob\n')
fprintf('%6.2f %0.2f\n',[theta0'; sum(w,1)])
fprintf('\n');
