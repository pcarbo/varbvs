%--------------------------------------------------------------------------
% varbvsprint.m: Print a summary of a model fitted by varbvs.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Print a summary of a Bayesian variable selection model fitted using 
%    variational approximation methods (see function varbvs).
%
% USAGE:
%    varbvsprint(fit)
%    varbvsprint(fit, c, nv, nr)
%
% INPUT ARGUMENTS:
% fit  Output of function varbvs.
% c    Compute c% credible intervals. By default, c = 0.95.
% nv   Show detailed statistics on nv variables. By default, nv = 5.
% nr   Draw nr samples for Monte Carlo estimates of the credible
%      intervals for the coefficients. Default is nr = 1000.
%
% DETAILS:
%    varbvsprint generates a four-part summary of the fitted Bayesian
%    variable selection model. The first part summarizes the analysis setup,
%    including the number of samples used to fit the model, the number of
%    covariates, and additional optimization settings. For the linear
%    regression model only, the "proportion of variance explained" (PVE) is
%    the posterior mean estimate of the proportion of variance in Y
%    explained by the variable selection model, along with a c% credible
%    interval of the PVE.
%
%    The second part summarizes the posterior distribution of the
%    hyperparameters by reporting the mean, the c% credible interval, and
%    the range of settings provided as input to varbvs (or NA when
%    options.sa or options.sigma is not set). These statistics are computed
%    by interpreting fit.logw as a vector of log-probabilities. (See
%    function varbvs for details, and how to adjust fit.logw to account for
%    different priors on the hyperparameters.) Summary statistics for
%    logodds are only shown if the prior log-odds is identical for all
%    variables.
%
%    The third part summarizes the variable selection results. Again, all
%    these statistics are computed by averaging over the hyperparameter
%    settings according to fit.logw. The number variables included at
%    different probability thresholds are given in the "count" row.
%
%    Finally, the fourth part gives more detailed statistics about the
%    variables that are most likely to be included in the regression model
%    (the 'prob.' column gives the posterior inclusion probability). Column
%    'coef' gives the posterior mean estimate of the regression coefficient
%    conditioned on being included in the model, and the next column gives
%    the c% credible interval. For the logistic regression model, the
%    coefficient is the same as the (natural) logarithm of the odds
%    ratio. For the linear regression model (family = 'gaussian'), the
%    posterior mean estimate of the proportion of variance in Y explained by
%    the variable ('PVE') is also provided.
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
%    varbvs, varbvsprint
%
% EXAMPLES:
%    See demo_qtl.m and demo_cc.m for examples.
%
function varbvsprint (fit, c, nv, nr)

  % Part of the varbvs package, https://github.com/pcarbo/varbvs
  %
  % Copyright (C) 2012-2017, Peter Carbonetto
  %
  % This program is free software: you can redistribute it under the
  % terms of the GNU General Public License; either version 3 of the
  % License, or (at your option) any later version.
  %
  % This program is distributed in the hope that it will be useful, but
  % WITHOUT ANY WARRANY; without even the implied warranty of
  % MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  % General Public License for more details.
  %

  % Get the number of variables (p) and number of candidate hyperparameter
  % settings (ns).
  p  = numel(fit.labels);
  ns = numel(fit.logw);

  % Take care of optional inputs.
  if nargin < 2
    c = 0.95;
  end

  % Show detailed statistics on nv variables. This cannot be greater than the
  % number of variables.
  if nargin < 3
    nv = 5;
  end
  nv = min(nv,p); 

  % Draw this many samples to compute the Monte Carlo estimates of the
  % credible intervals for the regression coefficients.
  if nargin < 4
    nr = 1000;
  end
  
  % (1) COMPUTE POSTERIOR STATISTICS
  % --------------------------------
  % Get the normalized (approximate) probabilities and the posterior
  % inclusion probabilities (PIPs).
  w   = fit.w;
  pip = fit.pip;
  
  % (2) SUMMARIZE ANALYSIS SETUP
  % ----------------------------
  fprintf('Summary of fitted Bayesian variable selection model:\n');
  fprintf('family:     %-8s',fit.family); 
  fprintf('   num. hyperparameter settings: %d\n',numel(fit.sa));
  fprintf('samples:    %-6d',fit.n); 
  fprintf('     iid variable selection prior: %s\n',tf2yn(fit.prior_same));
  fprintf('variables:  %-6d',p); 
  fprintf('     fit prior var. of coefs (sa): %s\n',tf2yn(fit.update_sa));
  fprintf('covariates: %-6d     ',size(fit.mu_cov,1));
  if strcmp(fit.family,'gaussian')
    fprintf('fit residual var. (sigma):    %s\n',tf2yn(fit.update_sigma));
  elseif strcmp(fit.family,'binomial')
    fprintf('fit approx. factors (eta):    %s\n',tf2yn(fit.optimize_eta));
  end
  fprintf('maximum log-likelihood lower bound: %0.4f\n',max(fit.logw));
  if strcmp(fit.family,'gaussian')
    x = sort(fit.model_pve);
    a = x(floor((0.5 - c/2)*length(x)));
    b = x(ceil((0.5 + c/2)*length(x)));
    fprintf('proportion of variance explained: ');
    fprintf('%0.3f [%0.3f,%0.3f]\n',mean(x),a,b);
  end

  % (3) SUMMARIZE RESULTS ON HYPERPARAMETERS
  % ----------------------------------------
  % Summarize the fitted residual variance parameter (sigma).
  fprintf('Hyperparameters: ');
  if ns == 1

    % Summarize the hyperparameter settings when there is only one
    % candidate setting.
    if strcmp(fit.family,'gaussian')
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
    if strcmp(fit.family,'gaussian')
      x = fit.sigma(:);
      if length(unique(x)) == 1
        fprintf('sigma   %8.3g NA                  %0.3g\n',x(1),x(1));
      else
        x0    = dot(w(:),x);
        [a b] = cred(x,x0,w,c);
        fprintf('sigma   %8.3g %-19s ',x0,sprintf('[%0.3g,%0.3g]',a,b));
        if fit.update_sigma
          fprintf('NA\n');
        else
          fprintf('%0.3g--%0.3g\n',min(x),max(x));
        end
      end
    end
 
    % Summarize the fitted prior variance parameter (sa).
    x = fit.sa(:);
    if length(unique(x)) == 1
      fprintf('sa      %8.3g NA                  %0.3g\n',x(1),x(1));
    else
      x0    = dot(w(:),x);
      [a b] = cred(x,x0,w,c);
      fprintf('sa      %8.3g %-19s ',x0,sprintf('[%0.3g,%0.3g]',a,b));
      if fit.update_sa
        fprintf('NA\n')
      else
        fprintf('%0.3g--%0.3g\n',min(x),max(x));
      end
    end
    
    % Summarize the fitted prior log-odds of inclusion (logodds).
    if (fit.prior_same)
      x = fit.logodds(:);
      if length(unique(x)) == 1
        fprintf('logodds %8.2f NA                  %0.2f\n',x(1),x(1));
      else
        x0    = dot(w(:),x);
        [a b] = cred(x,x0,w,c);
        fprintf('logodds %+8.2f %-19s (%+0.2f)--(%+0.2f)\n',x0,...
                sprintf('[%+0.2f,%+0.2f]',a,b),min(x),max(x));
      end
    end
  end
  
  % (4) SUMMARIZE VARIABLE SELECTION RESULTS
  % ----------------------------------------
  % Summarize the number of variables selected at different PIP thresholds.
  fprintf('Selected variables:\n');
  fprintf('prob. >0.10 >0.25 >0.50 >0.75 >0.90 >0.95\n');
  fprintf('count %5d %5d %5d %5d %5d %5d\n',sum(pip > 0.1),sum(pip > 0.25),...
          sum(pip > 0.5),sum(pip > 0.75),sum(pip > 0.9),sum(pip > 0.95));

  % Give more detailed statistics about the top nv variables by the
  % probability that they are included.
  [ans vars] = sort(-pip);
  vars       = vars(1:nv);
  vars       = vars(:)';
  fprintf('Top %d variables by inclusion probability:\n',nv);
  fprintf(' index variable   prob.');
  if strcmp(fit.family,'gaussian')
    fprintf('   PVE');
  end
  if strcmp(fit.family,'binomial')
    fprintf('   coef*   Pr(coef>%0.2f)\n',c);
  else
    fprintf('   coef    Pr(coef>%0.2f)\n',c);
  end
  for i = vars
    [a b] = varbvscoefcred(fit,i,c,nr);
    fprintf('%6d %-10s %0.3f',i,fit.labels{i},pip(i));
    if strcmp(fit.family,'gaussian')
      fprintf(' %0.3f',dot(w,fit.pve(i,:)));
    end
    fprintf(' %+7.3f [%+0.3f,%+0.3f]\n',fit.beta(i),a,b);
  end
  if strcmp(fit.family,'binomial')
    fprintf('*See "help varbvs" for interpreting coefficients in ');
    fprintf('logistic regression.\n');
  end

% ------------------------------------------------------------------
function y = tf2yn (x)
  if x
    y = 'yes';
  else
    y = 'no';
  end
