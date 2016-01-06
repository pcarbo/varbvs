% *** FIX THIS DESCRIPTION ***
%
% [ALPHA,MU1,MU2,XR] = VARBVSMIXUPDATE(X,SIGMA,SA1,SA2,LOGODDS,XY,D,...
% ALPHA0,MU10,MU20,XR0,I) runs a single iteration of the coordinate ascent
% updates to maximize the variational lower bound for Bayesian variable
% selection in linear regression. It adjusts the fully-factorized
% variational approximation to the posterior distribution of the
% coefficients in a linear regression model of a continuous outcome
% (quantitiative trait), with spike and slab priors on the coefficients.
%
% All inputs to this function are required. Input X is an N x P matrix of
% observations about the variables (or features), where N is the number of
% samples, and P is the number of variables. Input XY = X'*Y, where Y is the
% vector of observations about the outcome. Crucially, to account for an
% intercept, Y and X must be centered beforehand so that Y and each column
% of X has a mean of zero.
%
% This routine is implemented with the assumption that X is a single
% floating-point precision matrix (type HELP SINGLE), as opposed to the
% MATLAB default of double precision. This is useful for large data sets,
% because single precision requires half of the number of bits as double
% floating-point precision. If X is provided in another numerical
% representation, an error is reported.
%
% Inputs SIGMA, SA and LOGODDS specify the hyperparameters. SIGMA and SA are
% scalars. SIGMA specifies the variance of the residual, and SA*SIGMA is the
% prior variance of the regression coefficients. LOGODDS is the prior
% log-odds of inclusion for each variable. It is equal to LOGODDS =
% LOG(Q./(1-Q)), where Q is the prior probability that each variable is
% included in the linear model of Y. LOGODDS is a vector of length P.
%
% Inputs ALPHA0, MU0 are the current parameters of the variational
% approximation; under the variational approximation, the ith regression
% coefficient is normal with probability ALPHA0(i), and MU0(i) is the mean
% of the coefficient given that it is included in the model. Inputs XR0 and
% D must be XR0 = X*(ALPHA0.*MU0) and D = DIAG(X'*X).
%
% Input I specifies the order in which the coordinates are updated. It may
% be a vector of any length. Each entry of I must be an integer between 1
% and P.
%
% There are three outputs. Output vectors ALPHA and MU are the updated
% variational parameters, and XR = X*(ALPHA.*MU). The computational
% complexity of VARBVSUPDATE is O(N*LENGTH(I)).
function [alpha, mu1, mu2, Xr] = ...
        varbvsmixupdate (X, sigma, sa1, sa2, logodds, xy, d, ...
                         alpha0, mu10, mu20, Xr0, I)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % CHECK THE INPUTS.
  % X must be single precision.
  if ~isa(X,'double')
    error('Input argument X must be DOUBLE')
  end

  % Check inputs SIGMA, SA1 and SA2.
  if ~isscalar(sigma) | ~isscalar(sa1) | ~isscalar(sa2)
    error('Inputs SIGMA, SA and SA2 must be scalars');
  end

  % Check input LOGODDS.
  if isscalar(logodds)
    logodds = repmat(logodds,p,1);
  end
  if length(logodds) ~= p
    error('Input LOGODDS must be a scalar or a vector of length P');
  end
  
  % Check inputs XY and D.
  if length(xy) ~= p | length(d) ~= p
    error('Inputs XY and D must be vectors of length P');
  end

  % Check inputs ALPHA0, MU10 and MU20.
  if length(alpha0) ~= p | length(mu10) ~= p | length(mu20) ~= p
    error('Inputs ALPHA0, MU10 and MU20 must be vectors of length P');  
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
  [alpha mu1 mu2 Xr] = ...
      varbvsmixupdatematlab(X,double(sigma),double(sa1),double(sa2),...
                            double(logodds),double(xy),double(d),...
                            double(alpha0),double(mu10),double(mu20),...
                            double(Xr0),double(I-1));
