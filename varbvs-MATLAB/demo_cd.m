% This script fits the Bayesian variable selection model to identify genetic
% markers associated with Crohn's disease risk. The data consist of 442,001
% SNPs genotyped for 1,748 cases and 2,938 controls.
clear

% Initialize the random number generator. 
rng(1);

% LOAD GENOTYPE AND PHENOTYPE DATA
% --------------------------------
fprintf('LOADING DATA.\n');
load('/data/internal_restricted/carbonetto_2012_wtccc/MATLAB/cd.mat');
labels = strcat('rs',cellfun(@num2str,num2cell(labels),'UniformOutput',false));

% FIT VARIATIONAL APPROXIMATION TO POSTERIOR
% ------------------------------------------
% Fit the fully-factorized variational approximation to the posterior
% distribution of the coefficients for a logistic regression model of a
% binary outcome (case-control status), with spike and slab priors on the
% coefficients.
fprintf('FITTING MODEL TO DATA.\n')
fit = varbvs(X,[],y,labels,'binomial',struct('logodds',-6:0.25:-3));

% SUMMARIZE POSTERIOR DISTRIBUTION
% --------------------------------
fprintf('SUMMARIZING RESULTS.\n')
varbvsprint(fit,0.95,9);

% Compute "single-marker" posterior inclusion probabilities.
alpha = varbvsindep(fit,X,Z,y);

% Show two "genome-wide scans", one using the posterior inclusion
% probabilities (PIPs) computed in the joint analysis of all variables, and
% one using the PIPs that ignore correlations between the variables.
subplot(2,1,1);
varbvsplot(fit,struct('groups',chr,'n',9,'gap',5000));
subplot(2,1,2);
varbvsplot(fit,struct('groups',chr,'n',9,'gap',5000));
  
% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('varbvs_demo_cd.mat','fit','alpha','chr','pos','-v7.3');
