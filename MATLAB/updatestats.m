% [U,YHAT,XY,XU,D] = UPDATESTATS(X,Y,ETA) calculates useful quantities for
% updating the variational approximation to the logistic regression factors.
% Inputs X and Y specify the data (see VARBVSBIN), and ETA is the vector of
% free parameters. It is a column vector of length equal to the number of
% samples. This function should be called whenever the free parameters
% are modified.
%
% Outputs XY, XU and D are defined as XY = X'*YHAT, XU = X'*U and D =
% DIAG(X'*UHAT*X), where U = SLOPE(U). For a definition of vectors YHAT and
% matrix UHAT, see the Bayesian Analysis paper.
%
% The outputs may also be returned as a STRUCT.
function varargout = updatestats (X, y, eta)

  % Compute the slope of the conjugate.
  u = slope(eta);

  % Compute BETA0 and YHAT. See the journal paper for an explanation of
  % these two variables.
  beta0 = sum(y - 0.5)/sum(u);
  yhat  = y - 0.5 - beta0*u;

  % Here, I calculate XY = X'*YHAT as (YHAT'*X)' and XU = X'*U as (U'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(yhat'*X)';
  xu = double(u'*X)';

  % Compute the diagonal entries of X'*UHAT*X. For a definition of
  % UHAT, see the Bayesian Analysis journal paper.
  d = diagsq(X,u) - xu.^2/sum(u);

  % Take care of the variable output arguments.
  if nargout == 1
    varargout = { struct('u',u,'yhat',yhat,'xy',xy,'xu',xu,'d',d) };
  else
    varargout = { u yhat xy xu d };
  end