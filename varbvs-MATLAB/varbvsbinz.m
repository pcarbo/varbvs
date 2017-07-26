% [logw,sa,alpha,mu,s,eta] = varbvsbinz(X,Z,y,sa,logodds,...) implements
% the fully-factorized variational approximation for Bayesian variable
% selection in logistic regression, allowing for covariates. It is the same
% as varbvsbin, except that it allows for an additional set of covariates
% that are not subject to the same "spike-and-slab" priors as the other
% variables. The covariate data Z are specified as an n x m matrix, where n
% is the number of samples, and m is the number of covariates. This function
% is equivalent to varbvsbin when only one covariate is specified, the
% intercept, and Z = ones(n,1).
function [logw, err, sa, alpha, mu, s, eta] = ...
        varbvsbinz (X, Z, y, sa, logodds, alpha, mu, eta, tol, maxiter, ...
                   verbose, outer_iter, update_sa, optimize_eta, n0, sa0)
  
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

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % (1) INITIAL STEPS
  % -----------------
  % Input X must be single precision.
  if ~isa(X,'single')
    X = single(X);
  end

  % Compute a few useful quantities.
  Xr    = double(X*(alpha.*mu));
  stats = update_varbvsbinz_stats(X,Z,y,eta);
  s     = sa./(sa*stats.xdx + 1);

  % Initialize storage for outputs logw and err.
  logw = zeros(1,maxiter);
  err  = zeros(1,maxiter);

  % (2) MAIN LOOP
  % -------------
  % Repeat until convergence criterion is met, or until the maximum
  % number of iterations is reached.
  for iter = 1:maxiter
    
    % Save the current variational parameters and model parameters.
    alpha0 = alpha;
    mu0    = mu;
    s0     = s;
    eta0   = eta;
    sa0    = sa;

    % (2a) COMPUTE CURRENT VARIATIONAL LOWER BOUND
    % --------------------------------------------
    logw0 = int_logit(Z,y,stats,alpha,mu,s,Xr,eta) ...
            + int_gamma(logodds,alpha) ...
            + int_klbeta(alpha,mu,s,sa);

    % (2b) UPDATE VARIATIONAL APPROXIMATION
    % -------------------------------------
    % Run a forward or backward pass of the coordinate ascent updates.
    if mod(iter,2)
      i = 1:p;
    else
      i = p:-1:1;
    end
    [alpha mu Xr] = varbvsbinzupdate(X,sa,logodds,stats,alpha,mu,Xr,i);

    % (2c) UPDATE ETA
    % ---------------
    % Update the free parameters specifying the variational approximation
    % to the logistic regression factors.
    if optimize_eta
      eta   = update_eta(X,Z,y,betavar(alpha,mu,s),Xr,stats.d);
      stats = update_varbvsbinz_stats(X,Z,y,eta);
      s     = sa./(sa*stats.xdx + 1);
    end
   
    % (2d) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    % --------------------------------------------
    % Compute variational lower bound to marginal log-likelihood.
    logw(iter) = int_logit(Z,y,stats,alpha,mu,s,Xr,eta) ...
                 + int_gamma(logodds,alpha) ...
                 + int_klbeta(alpha,mu,s,sa);
    
    % (2e) UPDATE PRIOR VARIANCE OF REGRESSION COEFFICIENTS
    % -----------------------------------------------------
    % Compute the maximum a posteriori estimate of sa, if requested. Note
    % that we must also recalculate the variance of the regression
    % coefficients when this parameter is updated. 
    if update_sa
      sa = (sa0*n0 + dot(alpha,s + mu.^2))/(n0 + sum(alpha));
      s  = sa./(sa*stats.xdx + 1);
    end

    % (2f) CHECK CONVERGENCE
    % ----------------------
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small. If the variational bound decreases,
    % stop.
    err(iter) = max(abs(alpha - alpha0));
    if verbose
      if isempty(outer_iter)
        status = '';
      else
        status = sprintf('%05d ',outer_iter);
      end  
      status = [status sprintf('%05d %+13.6e %0.1e %06.1f      NA %0.1e',...
                               iter,logw(iter),err(iter),sum(alpha),sa)];
      fprintf(status);
      fprintf(repmat('\b',1,length(status)));
    end
    if logw(iter) < logw0
      logw(iter) = logw0;
      err(iter)  = 0; 
      sa         = sa0;
      alpha      = alpha0;
      mu         = mu0;
      s          = s0;
      eta        = eta0;
      break
    elseif err(iter) < tol
      break
    end
  end

  % Return the variational lower bound (logw) and "delta" in successive
  % iterates (err).
  logw = logw(1:iter);
  err  = err(1:iter);

% ----------------------------------------------------------------------
% diagprod(A,B) efficiently computes diag(A*B').
function y = diagprod (A, B)
  y = double(sum(A.*B,2));

% ----------------------------------------------------------------------
% Returns the M-step update for the parameters specifying the variational
% lower bound to the logistic regression factors, allowing for additional
% covariates.
function eta = update_eta (X, Z, y, v, Xr, d)

  % Compute muz, the posterior mean of the regression coefficients
  % corresponding to the covariates (u). Here, S is the posterior
  % covariance of u given beta.
  D   = diag(sparse(d));
  S   = inv(Z'*D*Z);
  muz = S*Z'*(y - 0.5 - d.*Xr);

  % Calculate the covariance between the coefficients u and beta.
  V = diag(sparse(v));
  C = -double((S*Z'*D)*X)*V;

  % This is the M-step update for the free parameters.
  U   = double((S*Z'*D)*X)';
  eta = sqrt((Z*muz + Xr).^2 + diagsq2(Z,S) + diagsq2(Z,U'*V*U) ...
             + diagsqt(X,v) + 2*diagprod(X*C',Z));

% -----------------------------------------------------------------------
function I = int_logit (Z, y, stats, alpha, mu, s, Xr, eta)

  % Get some of the statistics.
  yhat = stats.yhat;
  xdx  = stats.xdx;
  S    = stats.S;
  d    = stats.d;
  D    = diag(sparse(d));

  % Compute the variational approximation to the expectation of the
  % log-likelihood with respect to the approximate posterior distribution.
  I = sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 + logdet(S)/2 ...
      + qnorm(Z'*(y - 0.5),S)^2/2 + yhat'*Xr - qnorm(Xr,D)^2/2 ...
      + qnorm(Z'*D*Xr,S)^2/2 - xdx'*betavar(alpha,mu,s)/2;
