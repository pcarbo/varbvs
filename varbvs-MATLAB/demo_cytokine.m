% This script fits two variable selection models: the first ("null") model
% has a uniform prior for all variables (the 442,001 genetic markers), and
% the second model has higher prior probability for genetic markers near
% cytokine signaling genes. This main aim of this analysis is to assess
% support for enrichment of Crohn's disease risk factors near cytokine
% signaling genes; a large Bayes factor means greater support for the
% enrichment hypothesis. The data in this analysis consist of 442,001 SNPs
% genotyped for 1,748 cases and 2,938 controls.
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

% LOAD GENOTYPES, PHENOTYPES AND PATHWAY ANNOTATION
% -------------------------------------------------
fprintf('LOADING DATA.\n');
load('cytokine.mat');
load('cd.mat');
labels = strcat('rs',cellfun(@num2str,num2cell(labels),'UniformOutput',false));

% FIT VARIATIONAL APPROXIMATION
% -----------------------------
% Compute the variational approximation given the assumption that all
% variables (genetic markers) are, *a priori*, equally likely to be included
% in the model.
fprintf('FITTING NULL MODEL TO DATA.\n')
fit_null = varbvs(X,[],y,labels,'binomial',struct('logodds',-4));

% Compute the variational approximation given the assumption that genetic
% markers near cytokine signaling genes are more likely to be included in
% the model.
fprintf('FITTING PATHWAY ENRICHMENT MODEL TO DATA.\n')
logodds           = repmat(-4,442001,13);
logodds(a == 1,:) = repmat(-4 + (0:0.25:3),6711,1);
fit_cytokine = varbvs(X,[],y,labels,'binomial',...
                      struct('logodds',logodds,'alpha',fit_null.alpha,...
                             'mu',fit_null.mu,'eta',fit_null.eta,...
                             'optimize_eta',true));

% Compute the Bayes factor.
BF = bayesfactor(fit_null.logw,fit_cytokine.logw);

% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('~/data/varbvs_demo_cytokine.mat','fit_null','fit_cytokine','BF',...
     'a','chr','pos','-v7.3');

% Show two "genome-wide scans" from the multi-marker PIPs, with and without
% conditioning on enrichment of cytokine signaling genes.
subplot(2,1,1);
i = find(fit_null.alpha > 0.5 | fit_cytokine.pip > 0.5);
varbvsplot(fit_null,struct('groups',chr,'vars',i,'gap',5000));
ylabel('posterior probability');
subplot(2,1,2);
varbvsplot(fit_cytokine,struct('groups',chr,'vars',i,'gap',5000));
ylabel('posterior probability');
