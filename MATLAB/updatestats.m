% [U,YHAT,XY,XU,D] = UPDATESTATS(X,Y,ETA) calculates useful quantities for
% the posterior distribution of the additive genetic effects at each locus,
% under a variational approximation to the nonlinear terms in the logistic
% function. ETA is the vector of parameters for the variational lower bound
% to the logistic regression. It is a vector of length N.  This function
% should be called whenever the variational parameters (ETA) are updated.
%
% Input X is the genotype data. It is an N x P matrix, where N is the number
% of samples (individuals), and P is the number of variables ( SNPs). Y is
% the vector of quantitative trait data; it is a vector of length N. X and
% Y should not be centered.
%
% Outputs XY and XU are defined as XY = X'*YHAT and XU = X'*U. Output D is
% the diagonal of matrix UHAT; that is, D = diag(UHAT). For a definition of
% vectors U, YHAT and matrix UHAT, see the Bayesian Analysis journal paper.
%
% The outputs may also be returned as a STRUCT.
function varargout = updatestats (X, y, eta)

  % Compute the slope of the conjugate.
  u = slope(eta);

  % See the journal paper for a definition of BETA0 and YHAT.
  beta0 = sum(y - 0.5)/sum(u);
  yhat  = y - 0.5 - beta0*u;

  % Here, I calculate XY = X'*YHAT as (YHAT'*X)' and XU = X'*U as (U'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(yhat'*X)';
  xu = double(u'*X)';

  % Compute the diagonal entries of matrix UHAT. For a definition of
  % UHAT, see the journal paper.
  d = diagsq(X,u) - xu.^2/sum(u);

  % Take care of the variable output arguments.
  if nargout == 1
    varargout = { struct('u',u,'yhat',yhat,'xy',xy,'xu',xu,'d',d) };
  else
    varargout = { u yhat xy xu d };
  end