% [LOGW,ALPHA,MU,S,ETA] = MULTISNPBINZHYPER(X,Z,Y,H,THETA0,ALPHA,MU,ETA)
% runs the full variational inference procedure for Bayesian variable
% selection in logistic regression, allowing for covariates. It is the same
% as MULTISNPBINHYPER, except that it allows for an additional set of
% covariates that are not subject to the same "spike-and-slab" priors as the
% other variables. The covariate data Z are specified as an N x M matrix,
% where N is the number of samples, and M is the number of covariates. This
% function is equivalent to MULTISNPBINHYPER when only one covariate is
% specified, the intercept, and Z = ONES(N,1).
function [logw, alpha, mu, s, eta] = ...
        multisnpbinzhyper (X, Z, y, h, theta0, alpha, mu, eta)

  % Get the number of participants in the study (n), the number of SNPs
  % genotyped (p), and the number of combinations of the hyperparameters
  % (ns). 
  [n p] = size(X);
  ns    = numel(h);
  
  % First get the best initialization for the variational parameters.
  fprintf('Finding best initialization for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  [logw alpha mu s eta] = outerloophyperzbin(X,Z,y,alpha,mu,eta,h,theta0);

  % Choose an initialization common to all the runs of the coordinate ascent
  % algorithm. This is chosen from the hyperparameters with the highest
  % variational estimate of the posterior probability.
  [ans i] = max(logw(:));
  alpha   = repmat(alpha(:,i),1,ns);
  mu      = repmat(mu(:,i),1,ns);
  eta     = repmat(eta(:,i),1,ns);

  % Compute the unnormalized log-importance weights.
  fprintf('Computing importance weights for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  [logw alpha mu s eta] = outerloophyperzbin(X,Z,y,alpha,mu,eta,h,theta0);
