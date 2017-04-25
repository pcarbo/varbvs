% betavarmix(p,mu,s) returns variances of variables drawn from mixtures of
% normals. Each of the inputs is a n x k matrix, where n is the number of
% variables and k is the number of mixture components. Specifically,
% variable i is drawn from a mixture in which the jth mixture component is
% the univariate normal with mean mu(i,j) and variance s(i,j).
%
% Note that the following two lines should return the same result when k=2
% and the first component is the "spike" density with zero mean and
% variance.
%
%   y1 = betavar(p,mu,s)
%   y2 = betavarmix([1-p p],[zeros(n,1) mu],[zeros(n,1) s])
%
function y = betavarmix (p, mu, s)
  y = sum(p.*(s + mu.^2),2) - sum(p.*mu,2).^2;
