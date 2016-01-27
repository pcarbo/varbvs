% TO DO: Explain here what this script does.
clear

% SCRIPT PARAMETERS
% -----------------
logodds = (-6:0.25:-3);

% Initialize the random number generator. 
rng(1);

% LOAD GENOTYPE AND PHENOTYPE DATA
% --------------------------------
fprintf('LOADING DATA.\n');
load('/tmp/pcarbo/cd.mat');

% FIT VARIATIONAL APPROXIMATION TO POSTERIOR
% ------------------------------------------
% Fit the fully-factorized variational approximation to the posterior
% distribution of the coefficients for a logistic regression model of a
% binary outcome (case-control status), with spike and slab priors on the
% coefficients.
fprintf('FITTING MODEL TO DATA.\n')
fit = varbvs(X,[],y,labels,struct('logodds',logodds));
  
% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('/tmp/pcarbo/varbvs_demo_cd.mat','fit','-v7.3');
