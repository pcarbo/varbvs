%--------------------------------------------------------------------------
% varbvscoefcred.m: One-sentence summary of function goes here.
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
% OUTPUT ARGUMENTS:
% Description of output arguments goes here.
%
% DETAILS:
%    Detailed description of function goes here.
%
% LICENSE: GPL v3
%
% DATE: December 28, 2015
%
% AUTHORS:
%    List contributors here.
%
% REFERENCES:
%    List of references goes here.
%
% SEE ALSO:
%    List related functions here.
%
% EXAMPLES:
%    Give some examples here.
%
function [a, b] = varbvscoefcred (fit, vars, c, nr)

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

  % Compute the normalized (approximate) importance weights.
  w = normalizelogweights(fit.logw);

  % Initialize storage for the result.
  a = zeros(p,1);
  b = zeros(p,1);
  
  % Repeat for each selected variable.
  for i = 1:p
    j    = vars(i);
    k    = randtable(w,nr);
    x    = fit.mu(j,k) + sqrt(fit.s(j,k)) .* randn(1,nr);
    x    = sort(x);
    a(i) = x(floor((0.5 - c/2)*nr));
    b(i) = x(ceil((0.5 + c/2)*nr));
  end
  