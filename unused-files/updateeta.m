% UPDATEETA(X,Y,V,XR,D) returns the M-step update for the parameters
% specifying the variational lower bound to the logistic regression factors.
%
% The inputs are as follows. Input X is an N x P matrix of observations
% about the variables (or features), where N is the number of samples, and P
% is the number of variables. Y is the vector of observations about the
% binary trait; it is a vector of length N. Vector Y and the columns of
% matrix X must not be centered. Input XR must be equal to XR = X*R, where R
% is the posterior mean of the coefficients. Note that under the
% fully-factorized variational approximation, R = ALPHA.*MU; see VARBVSBIN
% for details.
%
% Input V is the posterior variance of the coefficients. For this update to
% be valid, it is required that the posterior covariance of the coefficients
% is equal to DIAG(V). Input D must be D = SLOPE(ETA), where ETA is the
% current estimate of the free parameters; see function SLOPE for details.
function eta = updateeta (X, y, v, Xr, d)
  
  % Compute MU0, the posterior mean of the intercept in the logistic
  % regression under the variational approximation. Here, A is the variance
  % of the intercept given the other coefficients.
  a   = 1/sum(d);
  mu0 = a*(sum(y - 0.5) - d'*Xr);

  % Compute S0, the (marginal) posterior variance of the intercept in the
  % logistic regression. Here, I calculate XD = X'*D as (D'*X)' to avoid
  % storing the transpose of X, since X may be large.
  xd = double(d'*X)';
  s0 = a*(1 + a*v'*(xd.^2));
  
  % Calculate the covariance between the intercept and coefficients.
  c = -a*xd.*v;
  
  % This is the M-step update for the free parameters.
  eta = sqrt((mu0 + Xr).^2 + s0 + diagsqt(X,v) + 2*double(X*c));
  