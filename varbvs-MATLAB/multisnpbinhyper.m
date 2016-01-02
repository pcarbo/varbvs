% [LOGW,ALPHA,MU,S,ETA] = MULTISNPBINHYPER(X,Y,H,THETA0) runs the full
% variational inference procedure for Bayesian variable selection in
% logistic regression (for a binary trait). The first part to the algorithm
% searches for a good initialization of the variational parameters.
% 
% This inference procedure involves an inner loop and an outer loop. The
% inner loop consists of running a coordinate ascent algorithm to tighten
% the variational lower bound given a setting of the hyperparameters. The
% outer loop computes importance weights for all combinations of the
% hyperparameters.
%
% Input X is the genotype data. It is an N x P matrix, where N is the number
% of samples (individuals), and P is the number of variables (genetic loci,
% or SNPs). Y is the vector of binary trait data, such as case-control
% status; it is a vector of length N. It is important to note that X and Y
% should *not* be centered. Instead, we account for the intercept while we
% update the variational approximation.
%
% Inputs H and THETA0 are the hyperparameter settings: H is the prior
% proportion of variance explained, and THETA0 is the (base 10) logarithm of
% the prior odds for inclusion; it is equal to THETA0 = LOG10(Q./(1-Q)),
% where Q is the prior probability that each SNP is included in the linear
% model of Y. These two inputs must be arrays of the same size. Note that a
% residual variance parameter SIGMA is not needed to model a binary
% trait. For each combination of the hyperparameters, we compute the
% unnormalized log-importance weight, and store the result in output LOGW.
%
% We assume a uniform prior on the proportion of variance explained, H. The
% prior on THETA0 is uniform over the interval given by the candidate
% values of THETA0. Assuming a uniform proposal distribution for the
% hyperparameters H and THETA0, the importance sampling procedure will
% recover the correct distribution of H and THETA0 as the number of Monte
% Carlo samples becomes large, and assuming the variational approximation
% closely matches the exact posterior. See PVE2SA for how to obtain the
% prior variance of the additive effects.
%
% Outputs ALPHA, MU and S are variational estimates of the posterior
% inclusion probabilities, and posterior means and variances of the additive
% effects (given that the variable is included in the model) for each
% setting of the hyperparameters. ALPHA, MU and S are each a matrix of
% dimension P x NS, where NS is the number of hyperparameter settings. Each
% column of ETA is the set of parameters specifying the variational
% approximation to logistic regression factors, conditioned on the
% hyperparameter settings. ETA is an N x NS matrix.
%
% [LOGW,ALPHA,MU,S,ETA] = MULTISNPBINHYPER(X,Y,H,THETA0,ALPHA0,MU0,ETA0)
% initializes the variational parameters for each combination of the
% hyperparameters, overriding a random initialization of these parameters.
function [logw, alpha, mu, s, eta] = ...
        multisnpbinhyper (X, y, h, theta0, alpha, mu, eta)

  % Get the number of participants in the study (n), the number of SNPs
  % genotyped (p), and the number of combinations of the hyperparameters
  % (ns). 
  [n p] = size(X);
  ns    = numel(h);

  % Set a random initialization of the variational parameters for each
  % combination of the hyperparameters, or use the initialization
  % provided by the input arguments.
  if exist('alpha','var')
    alpha = repmat(alpha,1,ns);
  else
    alpha = rand(p,ns);
    alpha = alpha ./ repmat(sum(alpha),p,1);
  end
  if exist('mu','var')
    mu = repmat(mu,1,ns);
  else
    mu = randn(p,ns);
  end
  if exist('eta','var')
    eta = repmat(eta,1,ns);
  else
    eta = ones(n,ns);
  end

  % First get the best initialization for the variational parameters.
  fprintf('Finding best initialization for %d combinations ',ns);
  fprintf('of hyperparameters.\n');
  [logw alpha mu s eta] = outerloophyperbin(X,y,alpha,mu,eta,h,theta0);

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
  [logw alpha mu s eta] = outerloophyperbin(X,y,alpha,mu,eta,h,theta0);
