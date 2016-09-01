%--------------------------------------------------------------------------
% varbvspredict.m: Make predictions from a model fitted by varbvs.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Predict outcomes (Y) given the observed variables (X) and observed
%    covariates (Z), and a model fitted by varbvs.
%
% USAGE:
%    varbvspredict(fit, X, Z)
%
% INPUT ARGUMENTS:
% fit  Output of function varbvs.
%
% X    n x p input matrix, in which p is the number of variables, and n
%      is the number of samples for which predictions will be made using
%      the fitted model. X cannot be sparse.
%
% Z    n x m covariate data matrix, where m is the number of covariates. Do
%      not supply an intercept as a covariate (i.e., a column of ones),
%      because an intercept is automatically included in the regression
%      model. For no covariates, set Z to the empty matrix [].
%
% OUTPUT: Vector containing the predicted outcomes for all samples. For
% family = 'binomial', all vector entries are 0 or 1.
%
% DETAILS:
%    For the logistic regression model, we do not provide classification
%    probabilities Pr(Y = 1 | X, Z) because these probabilities are not
%    necessarily calibrated under the variational approximation.
%   
%    The predictions are computed by averaging over the hyperparameter
%    settings, treating fit.logw as (unnormalized) log-marginal
%    probabilities. See varbvs for more details about correctly using
%    fit.logw for approximate numerical integration over the
%    hyperparameters, for example by treating these as importance
%    weights. 
%
% LICENSE: GPL v3
%
% DATE: February 19, 2016
%
% AUTHORS:
%    Algorithm was designed by Peter Carbonetto and Matthew Stephens.
%    R, MATLAB and C code was written by Peter Carbonetto.
%    Depts. of Statistics and Human Genetics, University of Chicago,
%    Chicago, IL, USA, and AncestryDNA, San Francisco, CA, USA
%
% REFERENCES:
%    P. Carbonetto, M. Stephens (2012). Scalable variational inference
%    for Bayesian variable selection in regression, and its accuracy in 
%    genetic association studies. Bayesian Analysis 7: 73-108.
%
% SEE ALSO:
%    varbvs
%
% EXAMPLES:
%    See demo_qtl.m and demo_cc.m for examples.
%
function y = varbvspredict (fit, X, Z)

  % Get the number of samples (n), variables (p) and hyperparameter
  % settings (ns).
  [n p] = size(X);
  ns    = numel(fit.logw);
    
  % Input X must be single precision, and cannot be sparse.
  if issparse(X)
    error('Input X cannot be sparse');
  end
  if ~isa(X,'single')
    X = single(X);
  end
  if (numel(fit.labels) ~= p)
    error('Inputs X and fit are not compatible');
  end

  % If input Z is not empty, it must be double precision, and must have as
  % many rows as X. Add an intercept to Z, and check the number of
  % covariates.
  if ~isempty(Z)
    if size(Z,1) ~= n
      error('Inputs X and Z do not match.');
    end
    Z = double(full(Z));
  end
  Z = [ones(n,1) Z];
  if (size(Z,2) ~= size(fit.mu_cov,1))
    error('Inputs Z and fit are not compatible');
  end

  % Compute the normalized (approximate) probabilities.
  w = normalizelogweights(fit.logw);

  % For each hyperparameter setting, and for each sample, compute the
  % posterior mean estimate of Y, and then average these estimates
  % over the hyperparameter settings. For the logistic regression, the
  % final "averaged" estimate is obtained by collecting the "votes"
  % from each hyperparameter setting, weighting the votes by the
  % marginal probabilities, and outputing the estimate that wins by
  % majority. The averaged estimate is computed this way because the
  % estimates conditioned on each hyperparameter setting are not
  % necessarily calibrated in the same way.
  Y = Z*fit.mu_cov + X*(fit.alpha.*fit.mu);
  if strcmp(fit.family,'gaussian')
    y = Y*w(:);
  elseif strcmp(fit.family,'binomial')
    y = round(round(sigmoid(Y))*w(:));
  else
    error('Invalid setting for fit.family');
  end