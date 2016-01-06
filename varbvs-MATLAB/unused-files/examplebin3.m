% In this smalls script, I test the "alternative" implementation of the
% variational approximation for logistic regression (VARBVSALTBIN) against
% the standard implementation (VARBVSBIN). Both procedures should give
% exactly the same solution.
clear

% SCRIPT PARAMETERS.
p  = 1000;  % Number of variables (SNPs).
n  = 2000;  % Number of samples (subjects).
na = 10;    % Number of SNPs that affect the outcome.
sa = 0.2;   % Standard deviation of log-odds ratios.
p1 = 0.4;   % Target proportion of subjects that are cases (y = 1).

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
beta        = sa * beta;
I           = find(beta ~= 0);
I           = I(randperm(na));
t           = na/2;
snps1       = I(1:t);
snps2       = I(t+1:na);
beta(snps2) = -beta(snps1);
maf(snps2)  = maf(snps1);
[X y]       = createbindata(maf,beta,logit(p1),n);

% COMPUTE VARIATIONAL ESTIMATES USING "varbvsbin".
fprintf('Computing variational estimates.\n');
[lnZ0 alpha0 mu0 s0 eta] = varbvsbin(X,y,sa^2,log(na/p));

% RECONSTRUCT VARIATIONAL ESTIMATES USING "varbvs".
options = struct('alpha',alpha0,'mu',mu0);
[lnZ alpha mu s] = varbvsaltbin(X,y,sa^2,log(na/p),eta,options);
