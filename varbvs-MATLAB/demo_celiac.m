% This script fits the Bayesian variable selection model to identify genetic
% markers associated with celiac disease risk. After removing all samples
% from the Finnish cohort, the data consist of 509,314 SNPs genotyped for
% 3,149 cases and 6,325 controls, and principal components (PCs) which are
% used as covariates in the logistic regression.
clear

% Initialize the random number generator. 
rng(1);

% LOAD GENOTYPE AND PHENOTYPE DATA
% --------------------------------
% Also load the principal components.
fprintf('LOADING DATA.\n');
load('/tmp/pcarbo/celiac_nomhc.mat');
load('/tmp/pcarbo/celiac_pc.mat');

% Select all samples *not* in the Finnish cohort.
i     = find(~strcmp(study,'Finnish'));
X     = X(i,:);
y     = y(i);
id    = id(i);
sex   = sex(i);
study = study(i);
pc    = pc(i,:);

% FIT VARIATIONAL APPROXIMATION TO POSTERIOR
% ------------------------------------------
% Fit the fully-factorized variational approximation to the posterior
% distribution of the coefficients for a logistic regression model of a
% binary outcome (case-control status), with spike and slab priors on the
% coefficients.
fprintf('FITTING MODEL TO DATA.\n')
fit = varbvs(X,pc(:,1:2),y,labels,'binomial',struct('logodds',-5.5:0.25:-3));
  
% SUMMARIZE POSTERIOR DISTRIBUTION
% --------------------------------
fprintf('SUMMARIZING RESULTS.\n')
% TO DO: Find an appropriate number of top-ranked variables to show in the
% summary.
varbvsprint(fit,0.95,10);

% TO DO: Compute "single-marker" posterior inclusion probabilities.

% TO DO: Show two "genome-wide scans", one using the multi-marker PIPs,
% and one using the single-marker PIPs. In the scan, label the top n SNPs
% by PIP.

% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('/tmp/pcarbo/varbvs_demo_celiac.mat','fit','chr','pos','-v7.3');
