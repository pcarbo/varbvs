% TO DO: Explain here what this script does.
clear

% Initialize the random number generator. 
rng(1);

% LOAD GENOTYPES, PHENOTYPES AND PATHWAY ANNOTATION
% -------------------------------------------------
fprintf('LOADING DATA.\n');
load('~/data/cytokine.mat');
load('/tmp/pcarbo/cd.mat');
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
logodds(a == 1,:) = repmat(0:0.25:3,6711,1);
fit_cytokine = varbvs(X,[],y,labels,'binomial',...
                      struct('logodds',logodds,'alpha',fit_null.alpha,...
                             'mu',fit_null.mu,'eta',fit_null.eta,...
                             'optimize_eta',true));

% Compute the Bayes factor.

% TO DO.

% SAVE RESULTS
% ------------
fprintf('SAVING RESULTS.\n');
save('/tmp/pcarbo/varbvs_demo_cytokine.mat','fit','a','BF','chr','pos',...
     '-v7.3');

% Show two "genome-wide scans" from the multi-marker PIPs, with and without
% conditioning on enrichment of cytokine signaling genes.
subplot(2,1,1);
w = normalizelogweights(fit_cytokine.logw);
i = find(fit_null.alpha > 0.5 | fit_cytokine.alpha*w(:) > 0.5);
varbvsplot(fit_null,struct('groups',chr,'vars',i,'gap',5000));
ylabel('posterior probability');
subplot(2,1,2);
varbvsplot(fit_cytokine,struct('groups',chr,'vars',i,'gap',5000));
ylabel('posterior probability');
