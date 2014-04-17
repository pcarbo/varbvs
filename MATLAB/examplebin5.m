% The goal of this script is to test a more efficient strategy for computing
% variational estimates of posterior probabilities and the variational lower
% bound when we only update the posterior probabilities for a small subset
% of SNPs.
clear

% SCRIPT PARAMETERS.
p  = 1000;  % Number of variables (SNPs).
n  = 4000;  % Number of samples (subjects).
na = 16;    % Number of SNPs that affect the outcome.
sa = 0.2;   % Standard deviation of log-odds ratios.
p1 = 0.4;   % Target proportion of subjects that are cases (y = 1).

% Regression coefficients corresponding to covariates (other than the
% intercept).
u = [-1 1]';

% Set the random number generator seed.
seed = 7;
rng(seed);

% GENERATE GENETIC EFFECTS.
% Create the additive effects on the outcome such that the first half of
% SNPs contains twice as many SNPs affecting the outcome than the rest of
% the SNPs.
t            = round(2/3*na);
[maf1 beta1] = createsnps(p/2,t);
[maf2 beta2] = createsnps(p/2,na-t);
maf          = [ maf1 maf2 ];
beta         = [ beta1; beta2 ];
clear maf1 maf2 beta1 beta2

% GENERATE DATA SET.
% To help ensure that the proportion of subjects with the disease ("cases")
% is equal to the target proportion (p1), I adjust the coefficients in the
% logistic regression so that every SNP that has an effect on disease status
% is matched with a SNP that has the same effect, but in the opposite
% direction. I also adjust the minor allele frequencies so that these
% "matched" SNPs have the same minor allele frequency.
fprintf('Creating data set.\n');
beta        = sa * beta;
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
Z     = [ones(n,1)] % X(:,p+(1:nz))];
X     = X(:,1:p);

% COMPUTE VARIATIONAL ESTIMATES USING UNIFORM PRIOR.
% First, compute the posterior inclusion probabilities with the same prior
% inclusion probability for all SNPs.
fprintf('Computing variational estimates using uniform prior.\n');
odds = repmat(2*na/(3*p),p,1);
[lnZ0 alpha0 mu0 s0 eta] = varbvszbin(X,Z,y,sa^2,log(odds));
fprintf('\n');

% COMPUTE VARIATIONAL ESTIMATES USING NON-UNIFORM PRIOR.
% Next, use a different prior inclusion probability for the first half of
% the SNPs (set "A"), and only update the posterior probabilities for the
% first half of SNPs.
fprintf('Computing variational estimates using non-uniform prior.\n');
t       = p/2;
A       = 1:t;
B       = t+1:p;
odds(A) = 2*odds(A);
options = struct('alpha',alpha0,'mu',mu0,'update_vars',1:t);
[lnZ alpha mu s] = varbvsaltzbin(X,Z,y,sa^2,log(odds),eta,options);
fprintf('\n');

% Finally, use a more efficient approach to update the variational
% estimates of the posterior probabilities for the first set of SNPs.
fprintf('Recomputing variational estimates for first half of SNPs.\n');
alpha2  = alpha0;
mu2     = mu0;
s2      = s0;
options = struct('alpha',alpha0(A),'mu',mu0(A),'Xb0',...
                 X(:,B)*(alpha0(B).*mu0(B)));
lnZa0   = varbvsaltzbin(X(:,A),Z,y,sa^2,log(2*na/(3*p)),eta,options);
fprintf('\n');
[lnZa1 alpha2(A) mu2(A) s2(A)] = ...
  varbvsaltzbin(X(:,A),Z,y,sa^2,log(4*na/(3*p)),eta,options);
