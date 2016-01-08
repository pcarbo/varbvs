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

  % Compute the posterior inclusion probabilities averaged over settings
  % of the hyperparameters.
  % TO DO.
  
  % Summarize the analysis setup.
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
  
  % Summarize the fitted residual variance parameter (sigma).
  fprintf('---\n');
  fprintf('         estimate Pr>%0.2f             candidates\n',c);
  if (fit.family == 'gaussian')
    x0    = dot(w(:),fit.sigma(:));
    [a b] = cred(fit.sigma,w,x0,c);
    fprintf('sigma:   %8.3g %-19s ',x0,sprintf('[%0.3g,%0.3g]',a,b));
    if fit.update_sigma
      fprintf('NA\n')
    else
      fprintf('%0.3g--%0.3g\n',min(fit.sigma(:)),max(fit.sigma(:)));
    end
  end

  % Summarize the fitted prior variance parameter (sa).
  x0    = dot(w(:),fit.sa(:));
  [a b] = cred(fit.sa,w,x0,c);
  fprintf('sa:      %8.3g %-19s ',x0,sprintf('[%0.3g,%0.3g]',a,b));
  if fit.update_sa
    fprintf('NA\n')
  else
    fprintf('%0.3g--%0.3g\n',min(fit.sa(:)),max(fit.sa(:)));
  end

  % Summarize the fitted prior log-odds of inclusion (logodds).
  if (fit.prior_same)
    x     = min(fit.logodds);
    x0    = dot(w(:),x);
    [a b] = cred(x,w,x0,c);
    fprintf('logodds: %+8.2f %-19s (%+0.2f)--(%+0.2f)\n',x0,...
            sprintf('[%+0.2f,%+0.2f]',a,b),min(x),max(x));
  end

  % 
  fprintf('---\n');

% Details to print:
%
%   - number of included variables at different probability thresholds.
%   - top n variables by PIP (+ mean and cred. int. of coefficients).
%   - Some stats on PVE---need to think about this.
%   - correlation between true and predicted Y.
%   

% ------------------------------------------------------------------
function y = tf2yn (x)
  if x
    y = 'yes';
  else
    y = 'no';
  end
