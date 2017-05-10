% [alpha,mu,Xr] = varbvsmixupdate(X,sigma,sa,q,xy,d,alpha0,mu0,Xr0,i)
% runs a single iteration of the coordinate ascent updates maximizing the
% variational lower bound for the linear regression model with a
% mixture-of-normals prior.
function [alpha, mu, Xr] = varbvsmixupdate (X, sigma, sa, q, xy, d, ...
                                            alpha0, mu0, Xr0, i)

  % I am adding this temporarily for testing.
  fast_version = false;
    
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
  if sum(i < 1 | i > p)
    error('Input i contains invalid variable indices');
  end

  % Initialize the outputs.
  alpha = alpha0;
  mu    = mu0;
  Xr    = Xr0;

  if (fast_version)
      
    % Execute the C routine. I subtract 1 from the indices because MATLAB 
    % arrays start at 1, but C arrays start at 0.
    [alpha mu Xr] = ...
      varbvsmixupdatemex(X,double(sigma),double(sa),double(q),double(xy),..
                         double(d),double(alpha0),double(mu0),double(Xr0),...
                         double(i-1));
  else
    i = i(:)';
    for j = i
      
      % Compute the variance of the regression coefficient conditioned on
      % being drawn from each of the mixture components. Note that the
      % variance corresponding to the first mixture component, the "spike",
      % should always be zero.
      s    = sigma*sa./(sa*d(j) + 1);
      s(1) = 0;
    
      % Update the variational estimate of the posterior mean for each
      % mixture component. Note that posterior mean corresponding to the
      % first component, the "spike", should always be zero.
      x       = double(X(:,j));
      r       = dot(alpha(j,:),mu(j,:));
      mu(j,:) = s/sigma * (xy(j) + d(j)*r - dot(x,Xr));
      mu(j,1) = 0;
    
      % Update the assignment probabilities for each of the mixture
      % components.
      SSR        = mu(j,:).^2./s;
      w          = log(q + eps) + (log(s./(sigma*sa)) + SSR)/2;
      w(1)       = log(q(1) + eps);
      alpha(j,:) = normalizelogweights(w);
    
      % Update Xr = X*r.
      rnew = dot(alpha(j,:),mu(j,:));
      Xr   = Xr + (rnew - r)*x;
    end
  end
