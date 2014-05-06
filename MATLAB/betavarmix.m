% DESCRIPTION GOES HERE.
function v = betavarmix (alpha, mu1, mu2, s1, s2)
  r = alpha.*mu1 + (1-alpha).*mu2;
  v = alpha.*(s1 + mu1.^2) + (1-alpha).*(s2 + mu2.^2) - r.^2;

