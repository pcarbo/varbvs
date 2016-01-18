% varbvspve(X,fit,nr) samples nr posterior estimates of the proportion of
% variance in Y explained by the Bayesian variable selection model fitted
% using a variational approximation. This function is only valid for the
% linear regression model (family = 'gaussian') with an intercept.
function pve = varbvspve (X, fit, nr)

  % Take care of the optional inputs.
  if nargin < 3
    nr = 1000;
  end

  % Get the number of variables.
  p = size(X,2);
    
  % Initialize storage for posterior estimates of the proportion of variance
  % explained.
  pve = zeros(nr,1);

  % Compute the normalized (approximate) importance weights.
  w = normalizelogweights(fit.logw);

  % For each sample, compute the proportion of variance explained.
  for i = 1:nr

    % Draw a hyperparameter setting from the posterior distribution.
    j = randtable(w);
    
    % Sample the region coefficients.
    b = fit.mu(:,j) + sqrt(fit.s(:,j)) .* randn(p,1);
    b = b .* (rand(p,1) < fit.alpha(:,j));

    % Compute the proportion of variance explained.
    sz     = var(X*b,1);
    pve(i) = sz/(sz + fit.sigma(j));
  end  
