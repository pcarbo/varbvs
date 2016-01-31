% TO DO: Add description of this function (see varbvscoefcred.m).
function [alpha, mu, s] = varbvsindep (fit, X, Z, y)

  % Get the number of samples (n), variables (p) and hyperparameter
  % settings (ns).
  [n p] = size(X);
  ns    = numel(fit.logw);

  % Input X must be single precision, and cannot be sparse. Note that
  % here I do not check the inputs as rigorously as function varbvs.
  if issparse(X)
    error('Input X cannot be sparse')
  end
  if ~isa(X,'single')
    X = single(X);
  end

  % If input Z is not empty, it must be double precision.
  if ~isempty(Z)
    Z = double(full(Z));
  end

  % Add intercept.
  Z = [ones(n,1) Z];

  % Input y must be a double-precision column vector.
  y = double(y(:));

  % If necessary, convert the prior log-odds to a p x ns matrix.
  if fit.prior_same
    logodds = repmat(logodds,p,1);
  end
  
  % Adjust the genotypes and phenotypes so that the linear effects of
  % the covariates are removed. This is equivalent to integrating out
  % the regression coefficients corresponding to the covariates with
  % respect to an improper, uniform prior; see Chipman, George and
  % McCulloch, "The Practical Implementation of Bayesian Model
  % Selection," 2001.
  if strcmp(family,'gaussian')
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
    if strcmp(family,'gaussian')
      [alpha(:,i) mu(:,i) s(:,i)] = ...
        varbvsnormindep(X,y,sigma(i),sa(i),log(10)*logodds(:,i));
    elseif strcmp(family,'binomial')
      [alpha(:,i) mu(:,i) s(:,i)] = ...
        varbvsbinzindep(X,Z,y,eta(:,i),sa(i),log(10)*logodds(:,i));
    end
  end
  