% [alpha,mu,Xr] = varbvsmixupdate(X,sigma,sa,q,xy,d,alpha0,mu0,Xr0,I)
% runs a single iteration of the coordinate ascent updates maximizing the
% variational lower bound for the linear regression model with a
% mixture-of-normals prior.
function [alpha, mu, Xr] = varbvsmixupdate (X, sigma, sa, q, xy, d, ...
                                            alpha0, mu0, Xr0, I)

  % Get the number of samples, (n) the number of variables (p), and the
  % number of mixture components including the "spike" (K)..
  [n p] = size(X);
  K     = length(q);

  % X should be single precision.
  if ~isa(X,'single')
    error('Input X should be SINGLE');
  end

  % Check input sigma.
  if ~isscalar(sigma)
    error('Input sigma should be a scalar');
  end

  % Check input xy and d.
  if ~(length(xy) == p & length(d) == p)
    error('Inputs xy and d should have length = size(X,2).');
  end
   
  % Check inputs alpha0 and mu0.
  if any([size(alpha0) size(mu0)] ~= [p K p K])
    error(cat(2,'Inputs alpha0 & mu0 should be p x K matrices, ',...
              'with p = size(X,2) and K = length(q).'));
  end
   
  % Check input Xr0.
  if length(Xr0) ~= n
    error('length(Xr0) should be equal to size(X,1)');
  end
   
  % Check input I.
  if sum(I < 1 | I > p)
    error('Input I contains invalid variable indices');
  end

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

