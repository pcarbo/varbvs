% BETAMEANMIX(P,MU,S) returns the mean of X, in which X is drawn from a
% mixture of two normals. Inputs MU1 and MU2 specify the means of the two
% normal densities, and P specifies the probability of drawing X from the
% first mixture component. Inputs P, MU1 and MU2 must be scalars, or arrays
% of the same dimension.
function r = betameanmix (p, mu1, mu2)
  r = p.*mu1 + (1-p).*mu2

