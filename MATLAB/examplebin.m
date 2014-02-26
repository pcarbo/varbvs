% In this small example, I simulate a data set from an idealized genetic
% association study in which the genetic markers (single nucleotide
% polymorphisms, or SNPs) are uncorrelated (i.e. no linkage disequilibrium),
% and I assess the performance of the variational approximation for mapping
% associations between the genetic markers and the binary trait (such as
% case-control status) in this setting.
%
% Note that I have tested this script in MATLAB 8.1, and it may not work
% correctly in earlier versions of MATLAB.
clear

% SCRIPT PARAMETERS.
p  = 1000;  % Number of variables (SNPs).
n  = 2e4;   % Number of samples (subjects).
na = 10;    % Number of SNPs that affect the outcome.
sb = 0.2;   % Standard deviation of log-odds ratios.
p1 = 0.4;   % Target proportion of subjects that are cases (y = 1).

% Candidate values for the prior proportion of variance explained (h), and
% prior log-odds for inclusion (theta0).
h      = (0.1:0.1:0.5)';
theta0 = (-2.5:0.25:-1.5)';

% Set the random number generator seed.
seed = 5;
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
[X y]       = createbindata(maf,beta,logit(p1),n);

% COMPUTE VARIATIONAL ESTIMATES.
fprintf('Computing variational estimates.\n');
[H THETA0] = ndgrid(h,theta0);
[logw alpha mu s eta] = multisnpbinhyper(X,y,H,THETA0);

% Compute the normalize importance weights.
w = normalizelogweights(logw);

% Calculate the posterior mean of the intercept.
mu0 = zeros(size(w));
r   = alpha .* mu;
for i = 1:numel(w)
  d      = slope(eta(:,i));
  mu0(i) = (sum(y - 0.5) - d'*X*r(:,i))/sum(d);
end
fprintf('Posterior mean of beta0: %0.2f\n',dot(mu0(:),w(:)));

% Calculate the posterior mean of the prior variance of the log-odds
% ratios (sa).
sx = sum(var1(X));
sa = pve2sa(sx,H,THETA0);
fprintf('Posterior mean of sa: %0.2f\n\n',dot(sa(:),w(:));

% Show the posterior distribution of the genome-wide log-odds (theta0). 
fprintf('Posterior of theta0:\n');
fprintf('theta0 prob\n')
fprintf('%6.2f %0.2f\n',[theta0'; sum(w,1)])
fprintf('\n');
