% *** DESCRIPTION OF FUNCTION GOES HERE ****
function varargout = updatestats2 (X, Z, y, eta)

  % Compute the slope of the conjugate.
  d = slope(eta);

  % Compute the posterior covariance of u (the regression coefficients
  % corresponding to the Z variables) given beta (the regression
  % coefficients corresponding to the X variables).
  D = diag(sparse(d));
  S = inv(Z'*D*Z);
  
  % Compute matrix DZR = D*Z*R, where R is an upper triangular matrix such
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