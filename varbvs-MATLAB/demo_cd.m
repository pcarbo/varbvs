% This script fits the Bayesian variable selection model to identify genetic
% markers associated with Crohn's disease risk. The data consist of 442,001
% SNPs genotyped for 1,748 cases and 2,938 controls.
clear

% Initialize the random number generator. 
rng(1);

% LOAD GENOTYPE AND PHENOTYPE DATA
% --------------------------------
fprintf('LOADING DATA.\n');
load('/tmp/pcarbo/cd.mat');
labels = strcat('rs',cellfun(@num2str,num2cell(labels),'UniformOutput',false));

% FIT VARIATIONAL APPROXIMATION TO POSTERIOR
% ------------------------------------------
% Fit the fully-factorized variational approximation to the posterior
% distribution of the coefficients for a logistic regression model of a
% binary outcome (case-control status), with spike and slab priors on the
% coefficients.
fprintf('FITTING MODEL TO DATA.\n')
fit = varbvs(X,[],y,labels,'binomial',struct('logodds',-6:0.25:-3));
  
% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('/tmp/pcarbo/varbvs_demo_cd.mat','fit','-v7.3');
