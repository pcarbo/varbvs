% [ALPHA,MU,XR] = VARBVSZBINUPDATE(X,SA,LOGODDS,STATS,ALPHA0,MU0,XR0,I) runs
% a single iteration of the coordinate ascent updates to maximize the
% variational lower bound for Bayesian variable selection in logistic
% regression, allowing for additional covariates. It is the same as
% VARBVSBINUPDATE, except that it allows for an arbitrary set of covariates.
%
% Input STATS is the STRUCT output from UPDATESTATS2.
function [alpha, mu, Xr] = varbvszbinupdate (X, sa, logodds, stats, ...
                                             alpha0, mu0, Xr0, I)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % CHECK THE INPUTS.
  % X must be single precision.
  if ~isa(X,'single')
    error('Input argument X must be SINGLE')
  end

  % Check input SA.
  if ~isscalar(sa)
    error('Input SA must be a scalar');
  end

  % Check input LOGODDS.
  if isscalar(logodds)
    logodds = repmat(logodds,p,1);
  end
  if length(logodds) ~= p
    error('Input LOGODDS must be a scalar or a vector of length P');
  end

  % Check input STATS.
  if ~isfield(stats,'d') | ~isfield(stats,'xy') | ~isfield(stats,'xdx') | ...
     ~isfield(stats,'dzr') 
    error('STATS must be a STRUCT with fields d, xy, xdx and dzr')
  end
  if length(stats.d) ~= n
    error('STATS.D must be a vector of length N');
  end
  if length(stats.xdx) ~= p | length(stats.xy) ~= p
    error('STATS.XDX and STATS.XY must be vectors of length P');
  end
  if size(stats.dzr,1) ~= n
    error('STATS.DZR must be a matrix with N rows');
  end
  stats.d   = double(stats.d);
  stats.xdx = double(stats.xdx);
  stats.xy  = double(stats.xy);
  stats.dzr = double(stats.dzr);

  % Check inputs ALPHA0 and MU0.
  if length(alpha0) ~= p | length(mu0) ~= p
    error('Inputs ALPHA0 and MU0 must be vectors of length P');  
  end

  % Check input XR0.
  if length(Xr0) ~= n
    error('Input XR0 must be a vector of length N');
  end

  % Check input I.
  if sum(I < 1 | I > p)
    error('Input I contains invalid variable indices');
  end

  % Execute the C routine. We need to subtract 1 from the indices because
  % MATLAB arrays start at one, but C arrays start at zero.
  [alpha mu Xr] = ...
      varbvszbinupdatematlab(X,double(sa),double(logodds),stats,...
                             double(alpha0),double(mu0),double(Xr0),...
                             double(I-1));
