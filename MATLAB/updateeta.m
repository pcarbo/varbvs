% *** FIX COMMENTS ***
%
% UPDATEETA(X,Y,V,XR,U) updates the parameters of the variational lower
% bound on the logistic regression model, assuming that the mean of the
% additive effects is R, and the covariance matrix of the additive effects
% is DIAG(V). Input XR must be equal to XR = X*R. Note that under our
% parameterization of the fully-factorized variational approximation, the
% mean of the additive effects is given by R = ALPHA.*MU; MU is the mean
% additive effect given that the SNP is included in the model.
%
% Input X is the genotype data. It is an N x P matrix, where N is the number
% of samples (individuals), and P is the number of variables (SNPs). Y is
% the vector of quantitative trait data; it is a vector of length N. It is
% important that X and Y not be centered. The return value is a vector of
% length N.
%
% For details on input U, see the help for functions SLOPE and UPDATESTATS.
function eta = updateeta (X, y, v, Xr, u)
  
  % Compute the posterior mean of the intercept in the logistic regression
  % (MU0) under the variational approximation. Here, A is the variance of
  % the intercept given the additive effects.
  a   = 1/sum(u);
  mu0 = a*(sum(y - 0.5) - u'*Xr);

  % Compute the (marginal) variance of the intercept in the logistic
  % regression. Here, I calculate XU = X'*U as (U'*X)' to avoid storing the
  % transpose of X, since X may be large.
  xu = double(u'*X)';
  s0 = a*(1 + a*v'*(xu.^2));
  
  % Calculate the covariance of the intercept and the additive effects.
  c = -a*xu.*v;
  
  % This is the M-step update for the variational parameters.
  eta = sqrt((mu0 + Xr).^2 + s0 + diagsq2(X,v) + 2*double(X*c));
  