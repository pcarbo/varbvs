% UPDATEETA(X,Z,Y,V,XR,D) returns the M-step update for the parameters
% specifying the variational lower bound to the logistic regression factors,
% allowing for additional covariates. It is the same as UPDATEETA, except
% for the inclusion of covariate data supplied in matrix Z.
function eta = updateetaz (X, Z, y, v, Xr, d)

  % Compute MUZ, the posterior mean of u, the regression coefficients
  % corresponding to the covariates. Here, S is the posterior covariance of
  % u given beta.
  D   = diag(sparse(d));
  S   = inv(Z'*D*Z);
  muz = S*Z'*(y - 0.5 - d.*Xr);

  % Calculate the covariance between the coefficients u and beta.
  V = diag(sparse(v));
  C = -double((S*Z'*D)*X)*V;

  % This is the M-step update for the free parameters.
  U   = double((S*Z'*D)*X)';
  eta = sqrt((Z*muz + Xr).^2 + diagsq2(Z,S) + diagsq2(Z,U'*V*U) ...
             + diagsqt(X,v) + 2*diagprod(X*C',Z));
