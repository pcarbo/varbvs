% TO DO: Add comments here explaining what function does.
function [alpha, mu, Xr] = varbvsmixupdate (X, sigma, sa, q, xy, d, ...
                                            alpha0, mu0, Xr0, I)

  % Get the number of samples (n) and variables (p).
  [n p] = size(X);

  % TO DO: Check inputs.

  % Initialize the outputs.
  alpha = alpha0;
  mu    = mu0;
  Xr    = Xr0;
  
  % Run the co-ordinate ascent updates.
  %
  % TO DO: Implement more efficient C routine.
  %
  I = I(:)';
  for i = I
      
    % Compute the variance of the regression coefficient conditioned on
    % being drawn from each of the mixture componentx. Note that the
    % variance corresponding to the first mixture component, the "spike",
    % should always be zero.
    s    = sigma*sa./(sa*d(i) + 1);
    s(1) = 0;
    
    % Update the variational estimate of the posterior mean for each
    % mixture component. Note that posterior mean corresponding to the
    % first component, the "spike", should always be zero.
    x       = double(X(:,i));
    r       = dot(alpha(i,:),mu(i,:));
    mu(i,:) = s/sigma * (xy(i) + d(i)*r - dot(x,Xr));
    mu(i,1) = 0;
    
    % Update the assignment probabilities for each of the mixture
    % components.
    SSR        = mu(i,:).^2./s;
    w          = log(q + eps) + (log(s./(sigma*sa)) + SSR)/2;
    w(1)       = log(q(1) + eps);
    alpha(i,:) = normalizelogweights(w);
    
    % Update Xr = X*r.
    rnew = dot(alpha(i,:),mu(i,:));
    Xr   = Xr + (rnew - r)*x;
  end

