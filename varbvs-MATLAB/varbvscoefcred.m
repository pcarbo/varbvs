%--------------------------------------------------------------------------
% varbvscoefcred.m: Compute credible intervals for posterior coefficients.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Compute Monte Carlo estimates of credible intervals for posterior mean
%    coefficients in the fitted variable selection model. This function is
%    used by varbvsprint to generate credible intervals for coefficients of
%    top-ranked variables.
%
% USAGE:
%    [a, b] = varbvscoefcred(fit, vars, c, nr)
%
% INPUT ARGUMENTS:
% fit   Output of function varbvs.
% vars  Compute c% credible intervals (CIs) for these variables.
%       By default, vars = 1:numel(fit.labels).
% c     Compute c% CIs. By default, c = 0.95.
% nr    Draw nr samples from posterior. Default is nr = 1000.
%
% OUTPUT ARGUMENTS:
% a     p x 1 array of CI lower bounds, where p = numel(vars).
% b     p x 1 array of CI upper bounds, where p = numel(vars). 
%
% DETAILS:
%    A c% credible interval (CI) is an interval [a,b] such that all values
%    between a and b account for >c% of the probability mass. This
%    description doesn't unique define the CI, and there are several ways to
%    arrive at a unique definition. Here, we define the CI by the (0.5 -
%    c/2)th and (0.5 + c/2)th quantiles of the posterior distribution.
%    Monte Carlo estimates of these quantiles are quickly computed by
%    drawing nr samples from the posterior.
%
%    Since these are Monte Carlo estimates, the output will differ
%    slightly for each call to varbvscoefcred. Larger nr will result in
%    more consistent outputs, but will also result in slower computation.
%
% LICENSE: GPL v3
%
% DATE: January 11, 2016
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
function [a, b] = varbvscoefcred (fit, vars, c, nr)

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

  % Get the number of hyperparameter settings.
  ns = numel(fit.logw);
  
  % Take care of optional inputs.
  if nargin < 2
    p    = length(fit.labels);
    vars = 1:p;
  else
    vars = vars(:)';
    p    = length(vars);
  end

  % By default, compute 95% credible intervals.
  if nargin < 3
    c = 0.95;
  end

  % By default, generate 1,000 random draws from the posterior
  % distribution. 
  if nargin < 4
    nr = 1000;
  end

  % Initialize storage for the result.
  a = zeros(p,1);
  b = zeros(p,1);
  
  % Repeat for each selected variable.
  for i = 1:p
    j    = vars(i);
    k    = randtable(fit.w,nr);
    x    = fit.mu(j,k) + sqrt(fit.s(j,k)) .* randn(1,nr);
    x    = sort(x);
    a(i) = x(floor((0.5 - c/2)*nr));
    b(i) = x(ceil((0.5 + c/2)*nr));
  end
  