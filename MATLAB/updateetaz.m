% *** DESCRIPTION OF FUNCTION GOES HERE ***
function eta = updateetaz (X, Z, y, v, Xr, d)

  % Compute MUZ, the posterior mean of u, the regression coefficients
  % corresponding to the covariates. Here, S is the posterior covariance of
  % u given beta.
  D   = diag(sparse(d));
  S   = inv(Z'*D*Z);
  muz = S*Z'*(y - 0.5 - d.*Xr);

  % Calculate the covariance between the coefficients u and beta.
  C = -double((S*Z'*D)*X)*diag(sparse(v));

  % This is the M-step update for the free parameters.
  % n   = length(d);
  % eta = zeros(n,1);
  % for i = 1:n
  %   x      = double(X(i,:)');
  %   z      = Z(i,:)';
  %   t      = double(X)'*D*(Z*(S*z));
  %   eta(i) = sqrt((z'*muz + Xr(i)).^2 + z'*S*z + t'*diag(sparse(v))*t ...
  %                 + dot(x.^2,v) + 2*z'*C*x);
  % end
  t   = double((S*Z'*D)*X)';
  eta = sqrt((Z*muz + Xr).^2 + diagsq2(Z,S) + diagsq2(Z,t'*diag(v)*t) ...
             + diagsqt(X,v) + 2*double(sum((Z*C).*X,2)));

