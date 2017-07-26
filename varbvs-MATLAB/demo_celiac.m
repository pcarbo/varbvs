% This script fits the Bayesian variable selection model to identify genetic
% markers associated with celiac disease risk. After removing all samples
% from the Finnish cohort, the data consist of 509,314 SNPs genotyped for
% 3,149 cases and 6,325 controls, and principal components (PCs) which are
% used as covariates in the logistic regression.
%
% Note that files celiac_nomhc.mat and celiac_pc.mat cannot be made publicly
% available due to data sharing restrictions, so this script is for viewing
% only.
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
% Also load the principal components.
fprintf('LOADING DATA.\n');
load('celiac_nomhc.mat');
load('celiac_pc.mat');

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
Z   = pc(:,1:2);
tic;
fit = varbvs(X,Z,y,labels,'binomial',struct('logodds',-5.5:0.25:-3));
r = toc;
fprintf('Model fitting took %0.2f minutes.\n',r/60);

% Compute "single-marker" posterior inclusion probabilities.
pip = varbvsindep(fit,X,Z,y) * fit.w(:);

% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('varbvs_demo_celiac.mat','fit','pip','chr','pos','r','-v7.3');

% SUMMARIZE POSTERIOR DISTRIBUTION
% --------------------------------
fprintf('SUMMARIZING RESULTS.\n')
varbvsprint(fit,0.95,14);

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
