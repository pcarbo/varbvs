% BETAVARMIX(P,MU1,MU2,S1,S2) returns the variance of X, in which X is drawn
% from a mixture of two normals. Inputs MU1, MU2, S1 and S2 specify the
% means and variances of the two normal densities, and P specifies the
% probability of drawing X from the first mixture component. Inputs P, MU1,
% MU2, S1 and S2 must be scalars, or arrays of the same dimension.
function v = betavarmix (p, mu1, mu2, s1, s2)
  v = p.*(s1 + mu1.^2) + (1-p).*(s2 + mu2.^2) - betameanmix(p,mu1,mu2).^2;

