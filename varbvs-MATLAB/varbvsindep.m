%--------------------------------------------------------------------------
% varbvsindep.m: Compute posterior statistics, ignoring correlations.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    This function computes the mean and variance of the coefficients, and
%    the posterior inclusion probabilities (PIPs), ignoring correlations
%    between variables. This is useful for inspecting or visualizing
%    groups of correlated variables (e.g., genetic markers in linkage
%    disequilibrium).
%
% USAGE:
%    [alpha, mu, s] = varbvsindep(fit, X, Z, y)
%
% INPUT ARGUMENTS:
% fit   Output of function varbvs.
%
% X     n x p input matrix, where n is the number of samples,
%       and p is the number of variables. X cannot be sparse.
%
% Z     n x m covariate data matrix, where m is the number of
%       covariates. Do not supply an intercept as a covariate
%       (i.e., a column of ones), because an intercept is
%       automatically included in the regression model. For no
%       covariates, set Z to the empty matrix [].
%
% y     Vector of length n containing observations of binary
%       (family = 'binomial') or continuous (family = 'gaussian')
%       outcome. For binary outcomes, all entries of y must be 
%       0 or 1.
%
% OUTPUT ARGUMENTS:
% alpha  Variational estimates of posterior inclusion probabilities.
% mu     Variational estimates of posterior mean coefficients.
% s      Variational estimates of posterior variances.
%
% DETAILS:
%    For the ith hyperparameter setting, alpha(:,i) is the variational
%    estimate of the posterior inclusion probability (PIP) for each
%    variable; mu(:,i) is the variational estimate of the posterior mean
%    coefficient given that it is included in the model; and s(:,i) is the
%    estimated posterior variance of the coefficient given that it is
%    included in the model.
%
% LICENSE: GPL v3
%
% DATE: January 31, 2016
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
function [alpha, mu, s] = varbvsindep (fit, X, Z, y)

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

  % If input Z is not empty, it must be double precision, and must have
  % as many rows as X.
  if ~isempty(Z)
    if size(Z,1) ~= n
      error('Inputs X and Z do not match.');
    end
    Z = double(full(Z));
  end

  % Add intercept.
  Z = [ones(n,1) Z];

  % Input y must be a double-precision column vector.
  y = double(y(:));

  % If necessary, convert the prior log-odds to a p x ns matrix.
  if fit.prior_same
    fit.logodds = repmat(fit.logodds,p,1);
  end
  
  % Adjust the genotypes and phenotypes so that the linear effects of
  % the covariates are removed. This is equivalent to integrating out
  % the regression coefficients corresponding to the covariates with
  % respect to an improper, uniform prior; see Chipman, George and
  % McCulloch, "The Practical Implementation of Bayesian Model
  % Selection," 2001.
  if strcmp(fit.family,'gaussian')
    if size(Z,2) == 1
      X = X - repmat(mean(X),length(y),1);
      y = y - mean(y);
    else

      % This should give the same result as centering the columns of X and
      % subtracting the mean from y when we have only one covariate, the
      % intercept.
      y = y - Z*((Z'*Z)\(Z'*y));
      X = X - Z*((Z'*Z)\(Z'*X));
    end
  end
  
  % Initialize storage for the outputs.
  alpha = zeros(p,ns);
  mu    = zeros(p,ns);
  s     = zeros(p,ns);

  % Calculate the mean (mu) and variance (s) of the coefficients given that
  % the coefficients are included in the model, and the posterior inclusion
  % probabilities (alpha), ignoring correlations between variables. Repeat
  % for each combination of the hyperparameters.
  for i = 1:ns
    if strcmp(fit.family,'gaussian')
      [alpha(:,i) mu(:,i) s(:,i)] = ...
        varbvsnormindep(X,y,fit.sigma(i),fit.sa(i),log(10)*fit.logodds(:,i));
    elseif strcmp(fit.family,'binomial')
      [alpha(:,i) mu(:,i) s(:,i)] = ...
        varbvsbinzindep(X,Z,y,fit.eta(:,i),fit.sa(i),log(10)*fit.logodds(:,i));
    end
  end
  