% TO DO: Explain here what this script does.
clear

% Initialize the random number generator. 
rng(1);

% LOAD GENOTYPE AND PHENOTYPE DATA
% --------------------------------
fprintf('LOADING DATA.\n');
load('/data/internal_restricted/carbonetto_2012_wtccc/MATLAB/cd.mat');
labels = strcat('rs',cellfun(@num2str,num2cell(labels),'UniformOutput',false));

% LOAD PATHWAY ANNOTATION
% -----------------------
load('~/data/cytokine.mat');

% FIT VARIATIONAL APPROXIMATION
% -----------------------------
% Compute the variational approximation given the assumption that all
% variables (genetic markers) are, *a priori*, equally likely to be included
% in the model.
fprintf('FITTING NULL MODEL TO DATA.\n')
% TO DO: FIX THIS.
fit_null = varbvs(X,[],y,labels,'binomial',struct('logodds',logodds));

% Compute the variational approximation given the assumption that genetic
% markers near cytokine signaling genes are more likely to be included in
% the model.
fprintf('FITTING PATHWAY ENRICHMENT MODEL TO DATA.\n')
% TO DO: FIX THIS.
fit_cytokine = varbvs(X,[],y,labels,'binomial',struct('logodds',logodds));

% TO DO: Compute the Bayes factor.

% TO DO: Show a "genome-wide scans" from the multi-marker PIPs, in which
% we condition on enrichment of cytokine signaling genes.

% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('varbvs_demo_cytokine.mat','fit','a','chr','pos','-v7.3');
