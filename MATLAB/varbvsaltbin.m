% [LNZ,ALPHA,MU,S] = VARBVSALTBIN(X,Y,SA,LOGODDS,ETA,OPTIONS) is exactly
% equivalent to VARBVSBIN with the same inputs, except that ETA, the vector
% of free parameters specifying the variational approximation to the
% logistic regression, must be chosen by hand. 
%
% Note that this implementation is potentially much slower than VARBVSBIN
% because it involves multiplying X times an N x N matrix, where N is the
% number of samples.
function [lnZ, alpha, mu, s] = varbvsaltbin (X, y, sa, logodds, eta, options)

  % Get the number of samples.
  n = length(y);

  % Take care of the optional inputs.
  if ~exist('options')
    options = [];
  end

  % Use the variational approximation to the logistic regression factors,
  % specified by the vector of parameters ETA, to transform the model to a
  % linear regression; that is, YHAT now acts as a quantitative trait.
  y    = y - 0.5;
  d    = slope(eta);
  D    = diag(d + 1e-6) - d*d'/sum(d);
  R    = chol(D);
  yhat = R'\(y - sum(y)/sum(d)*d);
  [lnZ alpha mu s] = varbvs(R*X,yhat,1,sa,logodds,options);

  % We need to modify the final expression for the variational lower
  % bound. 
  lnZ = lnZ + n/2*log(2*pi) - log(sum(d))/2 + yhat'*yhat/2 ...
      + sum(logsigmoid(eta)) + eta'*(d.*eta - 1)/2 ...
      + sum(y)^2/(2*sum(d));
