% TO DO: Add description of function here.
function y = varbvspredict (fit, X, Z)

  % Get the number of samples (n), variables (p) and hyperparameter
  % settings (ns).
  [n p] = size(X);
  ns    = numel(fit.logw);
    
  % Input X must be single precision, and cannot be sparse.
  if issparse(X)
    error('Input X cannot be sparse')
  end
  if ~isa(X,'single')
    X = single(X);
  end
  if (numel(fit.labels) ~= p)
    error('Inputs X and fit are not compatible')
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
    error('Inputs Z and fit are not compatible')
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