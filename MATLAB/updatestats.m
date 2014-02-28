% [D,YHAT,XY,XD,XDX] = UPDATESTATS(X,Y,ETA) calculates useful quantities for
% updating the variational approximation to the logistic regression factors.
% Inputs X and Y specify the data (see VARBVSBIN), and ETA is the vector of
% free parameters. It is a column vector of length equal to the number of
% samples. This function should be called whenever the free parameters
% are modified.
%
% Outputs XY, XD and XDX are defined as XY = X'*YHAT, XD = X'*D and XDX =
% DIAG(X'*DHAT*X), where D = SLOPE(ETA). For a definition of vectors YHAT
% and matrix DHAT, see the Bayesian Analysis paper.
%
% The outputs may also be returned as a STRUCT.
function varargout = updatestats (X, y, eta)

  % Compute the slope of the conjugate.
  d = slope(eta);

  % Compute BETA0 and YHAT. See the journal paper for an explanation of
  % these two variables.
  beta0 = sum(y - 0.5)/sum(d);
  yhat  = y - 0.5 - beta0*d;

  % Here, I calculate XY = X'*YHAT as (YHAT'*X)' and XD = X'*D as (D'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(yhat'*X)';
  xd = double(d'*X)';

  % Compute the diagonal entries of X'*DHAT*X. For a definition of
  % DHAT, see the Bayesian Analysis journal paper.
  xdx = diagsq(X,d) - xd.^2/sum(d);

  % Take care of the variable output arguments.
  if nargout == 1
    varargout = { struct('d',d,'yhat',yhat,'xy',xy,'xd',xd,'xdx',xdx) };
  else
    varargout = { d yhat xy xd xdx };
  end