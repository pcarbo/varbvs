% This script fits the Bayesian variable selection model to identify genetic
% markers associated with Crohn's disease risk. The data consist of 442,001
% SNPs genotyped for 1,748 cases and 2,938 controls.
%
% Note that file cd.mat cannot be made publicly available due to data
% sharing restrictions, so this script is for viewing only.
%
% Part of the varbvs package, https://github.com/pcarbo/varbvs
%
% Copyright (C) 2012-2017, Peter Carbonetto
%
% This program is free software: you can redistribute it under the
% terms of the GNU General Public License; either version 3 of the
% License, or (at your option) any later version.
%
% This program is distributed in the hope that it will be useful, but
% WITHOUT ANY WARRANY; without even the implied warranty of
% MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
% General Public License for more details.
%
clear

% Initialize the random number generator. 
rng(1);

% LOAD GENOTYPE AND PHENOTYPE DATA
% --------------------------------
fprintf('LOADING DATA.\n');
load('cd.mat');
labels = strcat('rs',cellfun(@num2str,num2cell(labels),'UniformOutput',false));

% FIT VARIATIONAL APPROXIMATION TO POSTERIOR
% ------------------------------------------
% Fit the fully-factorized variational approximation to the posterior
% distribution of the coefficients for a logistic regression model of a
% binary outcome (case-control status), with spike and slab priors on the
% coefficients.
fprintf('FITTING MODEL TO DATA.\n')
fit = varbvs(X,[],y,labels,'binomial',struct('logodds',-6:0.25:-3));

% Compute "single-marker" posterior inclusion probabilities.
pip = varbvsindep(fit,X,[],y) * fit.w(:);

% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('varbvs_demo_cd.mat','fit','pip','chr','pos','-v7.3');

% SUMMARIZE POSTERIOR DISTRIBUTION
% --------------------------------
fprintf('SUMMARIZING RESULTS.\n')
varbvsprint(fit,0.95,9);

% Show two "genome-wide scans", one using the posterior inclusion
% probabilities (PIPs) computed in the joint analysis of all variables, and
% one using the PIPs that ignore correlations between the variables. The
% latter is meant to look like a typical genome-wide "Manhattan" plot used
% to summarize the results of a genome-wide association study. Variables
% with PIP > 0.5 are highlighted.
i = find(fit.pip > 0.5);
subplot(2,1,1);
varbvsplot(fit,struct('groups',chr,'vars',i,'gap',5000));
ylabel('posterior probability');
subplot(2,1,2);
varbvsplot(fit,struct('groups',chr,'score',log10(pip + 0.001),'vars',i,...
                      'gap',5000));
ylabel('log10 posterior prob.');
