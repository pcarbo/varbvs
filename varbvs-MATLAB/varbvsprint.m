%--------------------------------------------------------------------------
% varbvsprint.m: Print a summary of a model fitted by varbvs.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Overview of function goes here.
%
% USAGE:
%    Summary of usage goes here.
%
% INPUT ARGUMENTS:
% Description of input arguments goes here.
%
% DETAILS:
%    Detailed description of function goes here.
%
% LICENSE: GPL v3
%
% DATE: January 10, 2016
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
%    varbvs.
%
% EXAMPLES:
%    See demo_qtl.m and demo_cc.m for examples.
%
function varbvsprint (fit, c, n, nr)

  % Get the number of variables (p) and number of candidate hyperparameter
  % settings (ns).
  p  = numel(fit.labels);
  ns = numel(fit.logw);

  % Take care of optional inputs.
  if nargin < 2
    c = 0.95;
  end

  % Show detailed statistics on n variables. Note that this cannot be
  % larger than the nubmer of variables.
  if nargin < 3
    n = 5;
  end
  n = min(n,p);

  % Draw this many samples to compute the Monte Carlo estimates of the
  % credible intervals for the regression coefficients.
  if nargin < 4
    nr = 1000;
  end
  
  % (1) COMPUTE POSTERIOR STATISTICS
  % --------------------------------
  % Compute the normalized (approximate) importance weights.
  w = normalizelogweights(fit.logw);

  % Compute the posterior inclusion probabilities (PIPs) and posterior mean
  % regression coefficients averaged over settings of the hyperparameters.
  PIP  = fit.alpha * w(:);
  beta = fit.mu    * w(:);
  
  % (2) SUMMARIZE ANALYSIS SETUP
  % ----------------------------
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
    fprintf('fit approx. factors (eta):    %s\n',tf2yn(fit.optimize_eta));
  end
  fprintf('maximum log-likelihood lower bound: %0.4f\n',max(fit.logw));
  if fit.family == 'gaussian'
    x  = sort(fit.model_pve);
    x0 = mean(x);
    a  = x(floor((0.5 - c/2)*length(x)));
    b  = x(ceil((0.5 + c/2)*length(x)));
    fprintf('proportion of variance explained: ');
    fprintf('%0.1f%% [%0.1f%%,%0.1f%%]\n',100*x0,100*a,100*b);
  end

  % (3) SUMMARIZE RESULTS ON HYPERPARAMETERS
  % ----------------------------------------
  % Summarize the fitted residual variance parameter (sigma).
  fprintf('Hyperparameters: ');
  if ns == 1

    % Summarize the hyperparameter settings when there is only one
    % candidate setting.
    if (fit.family == 'gaussian')
      fprintf('sigma=%0.3g ',fit.sigma);
    end
    fprintf('sa=%0.3g ',fit.sa);
    if (fit.prior_same)
      fprintf('logodds=%+0.2f',fit.logodds);
    end
    fprintf('\n');
  else
    fprintf('\n');
    fprintf('        estimate Pr>%0.2f             candidate values\n',c);
    if (fit.family == 'gaussian')
      x0    = dot(w(:),fit.sigma(:));
      [a b] = cred(fit.sigma,x0,w,c);
      fprintf('sigma   %8.3g %-19s ',x0,sprintf('[%0.3g,%0.3g]',a,b));
      if fit.update_sigma
        fprintf('NA\n')
      else
        fprintf('%0.3g--%0.3g\n',min(fit.sigma(:)),max(fit.sigma(:)));
      end
    end
 
    % Summarize the fitted prior variance parameter (sa).
    x0    = dot(w(:),fit.sa(:));
    [a b] = cred(fit.sa,x0,w,c);
    fprintf('sa      %8.3g %-19s ',x0,sprintf('[%0.3g,%0.3g]',a,b));
    if fit.update_sa
      fprintf('NA\n')
    else
      fprintf('%0.3g--%0.3g\n',min(fit.sa(:)),max(fit.sa(:)));
    end

    % Summarize the fitted prior log-odds of inclusion (logodds).
    if (fit.prior_same)
      x     = fit.logodds;
      x0    = dot(w(:),x);
      [a b] = cred(x,x0,w,c);
      fprintf('logodds %+8.2f %-19s (%+0.2f)--(%+0.2f)\n',x0,...
              sprintf('[%+0.2f,%+0.2f]',a,b),min(x),max(x));
    end
  end
  
  % (4) SUMMARIZE VARIABLE SELECTION RESULTS
  % ----------------------------------------
  % Summarize the number of variables selected at different PIP thresholds.
  fprintf('Selected variables:\n');
  fprintf('prob. >0.10 >0.25 >0.50 >0.75 >0.90 >0.95\n');
  fprintf('count %5d %5d %5d %5d %5d %5d\n',sum(PIP > 0.1),sum(PIP > 0.25),...
          sum(PIP > 0.5),sum(PIP > 0.75),sum(PIP > 0.9),sum(PIP > 0.95));

  % Give more detailed statistics about the top n variables by the
  % probability that they are included.
  [ans vars] = sort(-PIP);
  vars       = vars(1:n);
  vars       = vars(:)';
  fprintf('Top %d variables by inclusion probability:\n',n);
  fprintf('index variable   prob.');
  if fit.family == 'gaussian'
    fprintf(' -PVE-');
  end
  fprintf('   coef. Pr(coef.>%0.2f)\n',c);
  for i = vars
    [a b] = varbvscoefcred(fit,i,c,nr);
    fprintf('%5d %-10s %0.3f',i,fit.labels{i},PIP(i));
    if fit.family == 'gaussian'
      fprintf(' %04.1f%%',100*dot(w,fit.pve(i,:)));
    end
    fprintf(' %+7.3f [%+0.3f,%+0.3f]\n',beta(i),a,b);
  end

% ------------------------------------------------------------------
function y = tf2yn (x)
  if x
    y = 'yes';
  else
    y = 'no';
  end
