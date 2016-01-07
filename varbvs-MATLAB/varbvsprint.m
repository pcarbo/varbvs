% TO DO: Explain here what this function does.
function varbvsprint (fit, c)

  % Get the number of variables (p) and number of hyperparameter settings (ns).
  p  = numel(fit.labels);
  ns = numel(fit.logw);

  % Take care of optional inputs.
  if nargin < 2
    c = 0.95;
  end

  % Compute the normalized (approximate) importance weights.
  w = normalizelogweights(fit.logw);

  fprintf('Summary of fitted Bayesian variable selection model:\n')
  fprintf('family:     %-8s',fit.family); 
  fprintf('   num. hyperparameter settings: %d\n',numel(fit.sa));
  fprintf('samples:    %-6d',fit.num_samples); 
  fprintf('     iid variable selection prior: %s\n',tf2yn(fit.prior_same));
  fprintf('variables:  %-6d',p); 
  fprintf('     fit prior var. of coefs (sa): %s\n',tf2yn(fit.update_sa));
  fprintf('covariates: %-6d     ',fit.num_covariates);
  if fit.family == 'gaussian'
    fprintf('fit residual var. (sigma):    %s\n',tf2yn(fit.update_sigma));
  elseif fit.family == 'binomial'
    % TO DO: FIX THIS.
    fprintf('fit approx. factors (eta): %s\n',tf2yn(fit.optimize_eta));
  end
  fprintf('intercept:  %-3s        ',tf2yn(fit.intercept));
  fprintf('max. log-likelihood bound: %0.4f\n',max(fit.logw));
  fprintf('         estimate Pr>%0.2f\n',c);
  if (fit.family == 'gaussian')
    fprintf('sigma:   %0.2e\n',dot(w(:),fit.sigma(:)));
  end
    fprintf('sa:      %0.2e\n',dot(w(:),fit.sa(:)));
  if (fit.prior_same)
    fprintf('logodds: \n');
  end

% Details to print:
%
%   - range of candidate settings for sigma, sa and logodds
%   - sigma [+ credible interval]
%   - sa [+ credible interval]
%   - logodds [+ credible interval]
%   - number of included variables at different probability thresholds.
%   - top n variables by PIP (+ mean and cred. int. of coefficients).
%   - Some stats on PVE---need to think about this.
%   - correlation between true and predicted Y.
%   - max. logw
%   

% ------------------------------------------------------------------
function y = tf2yn (x)
  if x
    y = 'yes';
  else
    y = 'no';
  end
