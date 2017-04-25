% betavarmix(p,mu,s) returns variances of variables drawn from mixtures of
% normals. Each of the inputs is a p x k matrix, where p is the number of
% variables and k is the number of mixture components. Specifically,
% variable i is drawn from a mixture in which the jth mixture component is
% the univariate normal with mean mu(i,j) and variance s(i,j).
function y = betavarmix (p, mu, s)
  y = sum(p.*(s + mu.^2),2) - sum(p.*mu,2).^2;
