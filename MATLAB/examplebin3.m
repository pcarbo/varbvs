% PROVDE SHORT DESCRIPTION OF MATLAB SCRIPT HERE.
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
d    = slope(eta);
D    = diag(d) - d*d'/sum(d);
R    = chol(D);
Xhat = R*X;
yhat = y - 0.5;
yhat = R'\(yhat - sum(yhat)/sum(d)*d);
options = struct('alpha',alpha0,'mu',mu0);
[lnZ alpha mu s] = varbvs(Xhat,yhat,1,sa^2,log(na/p),options);
lnZ = lnZ + n/2*log(2*pi) - log(sum(d))/2 + yhat'*yhat/2 ...
      + sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 ...
      + sum(y - 0.5)^2/(2*sum(d));