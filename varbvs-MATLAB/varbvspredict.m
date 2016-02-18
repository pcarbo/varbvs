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

  % TO DO: Add check for number of covariates.

  % Compute the normalized (approximate) probabilities.
  w = normalizelogweights(fit.logw);

  % For each hyperparameter setting, and for each sample, compute the
  % posterior mean estimate of Y, then return the final estimate averaged
  % over the hyperparameter settings.
  y = zeros(n,ns);
  for i = 1:ns
    y(:,i) = Z*fit.muz(:,i) + X*(fit.alpha(:,i).*fit.mu(:,i));
  end
  if strcmp(fit.family,'gaussian')
    y = y * w(:);
  elseif strcmp(fit.family,'binomial')
    y = round(sigmoid(y) * w(:));
  end
