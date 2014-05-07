% DESCRIPTION OF FUNCTION GOES HERE.
function [lnZ, alpha, mu1, mu2, s1, s2] = ...
    varbvsmix (X, y, sigma, sa1, sa2, logodds)

  % Convergence is reached when the maximum relative distance between
  % successive updates of the variational parameters is less than this
  % quantity.
  tolerance = 1e-4;  

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % CHECK INPUTS.
  % TO DO.

  % TAKE CARE OF OPTIONAL INPUTS.
  % TO DO.

  % Set initial estimates of variational parameters.
  alpha = rand(p,1);
  alpha = alpha / sum(alpha);
  mu1   = sqrt(sigma*sa1) * randn(p,1);
  mu2   = sqrt(sigma*sa2) * randn(p,1);

  % INITIAL STEPS.
  % Compute a few useful quantities. Here I calculate X'*Y as (Y'*X)' to
  % avoid storing the transpose of X, since X may be large.
  xy = double(y'*X)';
  d  = diagsq(X);
  Xr = double(X*(alpha.*mu1 + (1-alpha).*mu2));
  
  % Calculate the variances of the coefficients.
  s1 = sa1*sigma./(sa1*d + 1);
  s2 = sa2*sigma./(sa2*d + 1);

  % MAIN LOOP.
  % Repeat until convergence criterion is met.
  lnZ  = -Inf;
  iter = 0;
  fprintf('       variational    max. incl max.      \n');
  fprintf('iter   lower bound  change vars E[b] sigma\n');
  while true

    % Go to the next iteration.
    iter = iter + 1;

    % Save the current variational parameters and lower bound.
    alpha0  = alpha;
    mu10    = mu1;
    mu20    = mu2;
    lnZ0    = lnZ;
    params0 = [ alpha; alpha.*mu1; (1-alpha).*mu2 ];

    % UPDATE VARIATIONAL APPROXIMATION.
    % Run a forward or backward pass of the coordinate ascent updates.
    if isodd(iter)
      I = 1:p;
    else
      I = p:-1:1;
    end
    for i = 1:p
  
      % Update the variational estimates of the posterior means.
      r      = alpha(i)*mu1(i) + (1-alpha(i))*mu2(i);
      mu1(i) = s1(i)/sigma * (xy(i) + d(i)*r - dot(X(:,i),Xr));
      mu2(i) = s2(i)/sigma * (xy(i) + d(i)*r - dot(X(:,i),Xr));
  
      % Update the variational estimate of the posterior inclusion
      % probability.
      SSR1     = mu1(i)^2/s1(i);
      SSR2     = mu2(i)^2/s2(i);
      alpha(i) = sigmoid(logodds + (log(s1(i)*sa2/(s2(i)*sa1)) + SSR1-SSR2)/2);
  
      % Update Xr = X*r.
      rnew = alpha(i)*mu1(i) + (1-alpha(i))*mu2(i);
      Xr   = Xr + (rnew - r)*X(:,i);
    end

    % UPDATE RESIDUAL VARIANCE.
    sigma = (norm(y - Xr)^2 + d'*betavarmix(alpha,mu1,mu2,s1,s2) ...
             + alpha'*(s1 + mu1.^2)/sa1 ...
             + (1-alpha)'*(s2 + mu2.^2)/sa2)/(n + p);
    s1 = sa1*sigma./(sa1*d + 1);
    s2 = sa2*sigma./(sa2*d + 1);

    % COMPUTE VARIATIONAL LOWER BOUND.
    % Compute the lower bound to the marginal log-likelihood.
    lnZ = - n/2*log(2*pi*sigma) - norm(y - Xr)^2/(2*sigma) ...
          - d'*betavarmix(alpha,mu1,mu2,s1,s2)/(2*sigma) ...
          + intgamma(logodds,alpha) ...
	  + intklbeta(alpha,mu1,s1,sigma*sa1) ...
          + intklbeta(1-alpha,mu2,s2,sigma*sa2) ...
          + alpha'*log(alpha + eps) + (1 - alpha)'*log(1 - alpha + eps);

    % CHECK CONVERGENCE.
    % Print the status of the algorithm and check the convergence criterion.
    % Convergence is reached when the maximum relative difference between
    % the parameters at two successive iterations is less than the specified
    % tolerance, or when the variational lower bound has decreased. I ignore
    % parameters that are very small.
    params = [ alpha; alpha.*mu1; (1-alpha).*mu2 ];
    I      = find(abs(params) > 1e-6);
    err    = relerr(params(I),params0(I));
    fprintf('%4d %+13.6e %0.1e %4d %0.2f %5.2f\n',iter,lnZ,max(err),...
            round(sum(alpha)),max(abs(alpha.*mu1)),sqrt(sigma));
    if lnZ < lnZ0
      alpha = alpha0;
      mu1   = mu10;
      mu2   = mu20;
      lnZ   = lnZ0;
      break
    elseif max(err) < tolerance
      break
    end
  end