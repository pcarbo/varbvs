% [logw,sa,alpha,mu,s,eta] = varbvsbin(X,y,sa,logodds,...) implements the
% fully-factorized variational approximation for Bayesian variable selection
% in logistic regression. It finds the "best" fully-factorized variational
% approximation to the posterior distribution of the coefficients in a
% logistic regression model of a binary outcome or trait, with spike and
% slab priors on the coefficients. By "best", we mean the approximating
% distribution that locally minimizes the K-L divergence between the
% approximating distribution and the exact posterior.
%
% Input X is an n x p matrix variable (or feature) observations, where n is
% the number of samples, and p is the number of variables. Input y contains
% samples of the binary outcome; it is a vector of length n.
%
% Inputs sa and logodds are the hyperparameters. Scalar sa is the prior
% variance of the coefficients. Input logodds is the prior log-odds of
% inclusion for each variable. Note that the prior log-odds here is defined
% with respect to the *natural* logarithm, whereas in function varbvs the
% prior log-odds is defined with respect to the base-10 logarithm, so a
% scaling factor of log(10) is needed to convert from the latter to the
% former. Also, note that the residual variance parameter (sigma) is not
% needed to model a binary outcome.
%
% Output logw is the variational estimate of the marginal log-likelihood
% given the hyperparameters at each iteration of the co-ordinate ascent
% optimization procedure. Output err is the maximum difference between the
% approximate posterior probabilities (alpha) at successive iterations.
% Outputs alpha, mu and s are the parameters of the variational
% approximation and, equivalently, variational estimates of posterior
% quantites: under the variational approximation, the ith regression
% coefficient is normal with probability alpha(i); mu(i) and s(i) are the
% mean and variance of the coefficient given that it is included in the
% model. Output eta is the vector of free parameters that specify the
% variational approximation to the likelihood factors in the logistic
% regression. Finally, output mu0 is the posterior mean of the intercept.
function [logw, err, sa, alpha, mu, s, eta, mu0] = ...
        varbvsbin (X, y, sa, logodds, alpha, mu, eta, tol, maxiter, ...
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
  stats = update_varbvsbin_stats(X,y,eta);
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
    % Compute variational lower bound to marginal log-likelihood.
    logw0 = int_logit(y,stats,alpha,mu,s,Xr,eta) ...
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
    [alpha mu Xr] = varbvsbinupdate(X,sa,logodds,stats,alpha,mu,Xr,i);

    % (2c) UPDATE ETA
    % ---------------
    % Update the free parameters specifying the variational approximation
    % to the logistic regression factors.
    if optimize_eta
      eta   = update_eta(X,y,betavar(alpha,mu,s),Xr,stats.d);
      stats = update_varbvsbin_stats(X,y,eta);
      s     = sa./(sa*stats.xdx + 1);
    end

    % (2d) COMPUTE UPDATED VARIATIONAL LOWER BOUND
    % --------------------------------------------
    % Compute variational lower bound to marginal log-likelihood.
    logw(iter) = int_logit(y,stats,alpha,mu,s,Xr,eta) ...
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
% update_eta(X,y,v,Xr,d) returns the M-step update for the parameters
% specifying the variational lower bound to the logistic regression factors.
% Input Xr must be Xr = X*r, where r is the posterior mean of the
% coefficients. Note that under the fully-factorized variational
% approximation, r = alpha.*mu. Input v is the posterior variance of the
% coefficients. For this update to be valid, it is required that the
% posterior covariance of the coefficients is equal to diag(v). Input d must
% be d = slope(eta); see function slope for details.
function eta = update_eta (X, y, v, Xr, d)
  
  % Compute mu0, the posterior mean of the intercept in the logistic
  % regression under the variational approximation. Here, a is the
  % conditional variance of the intercept given the other coefficients.
  a   = 1/sum(d);
  mu0 = a*(sum(y - 0.5) - d'*Xr);

  % Compute s0, the (marginal) posterior variance of the intercept in the
  % logistic regression. Here, I calculate xd = X'*d as (d'*X)' to avoid
  % storing the transpose of X, since X may be large.
  xd = double(d'*X)';
  s0 = a*(1 + a*v'*(xd.^2));
  
  % Calculate the covariance between the intercept and coefficients.
  c = -a*xd.*v;
  
  % This is the M-step update for the free parameters.
  eta = sqrt((mu0 + Xr).^2 + s0 + diagsqt(X,v) + 2*double(X*c));

% ----------------------------------------------------------------------
% int_logit(y,stats,alpha,mu,s,Xr,eta) computes an integral that appears in
% the variational lower bound of the marginal log-likelihood for the
% logistic regression model. This integral is an approximation to the
% expectation of the logistic regression log-likelihood taken with respect
% to the variational approximation. Input Xr must be equal to Xr = X*r,
% where r is the posterior mean of the coefficients. Note that under the
% fully-factorized variational approximation, r = alpha.*mu.
function I = int_logit (y, stats, alpha, mu, s, Xr, eta)

  % Get some of the statistics.
  yhat = stats.yhat;
  xdx  = stats.xdx;
  d    = stats.d;
  D    = diag(sparse(d));

  % Get the variance of the intercept given the other coefficients.
  a = 1/sum(d);

  % Compute the variational approximation to the expectation of the
  % log-likelihood with respect to the variational approximation.
  I = sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 + log(a)/2 ...
      + a*sum(y - 0.5)^2/2 + yhat'*Xr - qnorm(Xr,D)^2/2 ...
      + a*(d'*Xr)^2/2 - xdx'*betavar(alpha,mu,s)/2;
