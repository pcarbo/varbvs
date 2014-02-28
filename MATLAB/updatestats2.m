% [S,D,YHAT,XY,XD,XDX,DZR] = UPDATESTATS2(X,Z,Y,ETA) calculates useful
% quantities for updating the variational approximation to the logistic
% regression factors, allowing for covariates. Inputs X, Z and Y specify the
% data (see VARBVSZBIN), and ETA is the vector of free parameters. It is a
% column vector of length equal to the number of samples. This function
% should be called whenever the free parameters are modified.
%
% The outputs are defined as XY = X'*YHAT, XD = X'*D, XDX = DIAG(X'*DHAT*X),
% DZR = D*Z*R' and S = INV(Z'*D*Z) is the posterior covariance of the
% covariate effects given the other effects, where D = SLOPE(ETA), and R is
% the upper triangular matrix such that R'*R = S. For a definition of
% vectors YHAT and matrix DHAT, see the Bayesian Analysis paper.
%
% The outputs may also be returned as a STRUCT.
function varargout = updatestats2 (X, Z, y, eta)

  % Compute the slope of the conjugate.
  d = slope(eta);

  % Compute the posterior covariance of u (the regression coefficients
  % corresponding to the Z variables) given beta (the regression
  % coefficients corresponding to the X variables).
  D = diag(sparse(d));
  S = inv(Z'*D*Z);
  
  % Compute matrix DZR = D*Z*R', where R is an upper triangular matrix such
  % that R'*R = S.
  R   = chol(S);
  dzr = D*(Z*R');

  % Compute YHAT. 
  yhat = y - 0.5 - dzr*R*(Z'*(y - 0.5));

  % Here, I calculate XY = X'*YHAT as (YHAT'*X)' and XD = X'*D as (D'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(yhat'*X)';
  xd = double(d'*X)';

  % Compute the diagonal entries of X'*DHAT*X. For a definition of DHAT, see
  % the Bayesian Analysis journal paper.
  xdx = diagsq(X,d) - diagsq(dzr'*X);

  % Take care of the variable output arguments.
  if nargout == 1
    varargout = { struct('S',S,'d',d,'yhat',yhat,'xy',xy,'xd',xd,...
                         'xdx',xdx,'dzr',dzr) };
  else
    varargout = { S d yhat xy xd xdx dzr };
  end