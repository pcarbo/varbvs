% [LNZ,ALPHA,MU,S] = VARBVSALTZBIN(X,Z,Y,SA,LOGODDS,ETA,OPTIONS) is exactly
% equivalent to VARBVSZBIN with the same inputs, except that ETA, the vector
% of free parameters specifying the variational approximation to the
% logistic regression, must be chosen by hand. 
%
% Note that this implementation is potentially much slower than VARBVSZBIN
% because it involves multiplying X times an N x N matrix, where N is the
% number of samples.
function [lnZ, alpha, mu, s] = varbvsaltzbin (X, Z, y, sa, logodds, ...
                                              eta, options)

  % Get the number of samples.
  n = length(y);

  % Take care of the optional inputs.
  if ~exist('options')
    options = [];
  end

  % Compute the posterior covariance of u (the regression coefficients
  % corresponding to the Z variables) given beta (the regression
  % coefficients corresponding to the X variables).
  d  = slope(eta);
  DZ = diag(sparse(d))*Z;
  S  = inv(Z'*DZ);

  % Use the variational approximation to the logistic regression factors,
  % specified by the vector of parameters ETA, to transform the model to a
  % linear regression; that is, YHAT now acts as a quantitative trait.
  y    = y - 0.5;
  D    = diag(d + 1e-6) - DZ*S*DZ';
  R    = chol(D);
  yhat = R'\(y - DZ*S*(Z'*y));
  [lnZ alpha mu s] = varbvs(R*X,yhat,1,sa,logodds,options);

  % We need to modify the final expression for the variational lower
  % bound. 
  lnZ = lnZ + n/2*log(2*pi) + logdet(S)/2 + yhat'*yhat/2 ... 
        + sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 ...
        + qnorm(Z'*y,S)^2/2;
  