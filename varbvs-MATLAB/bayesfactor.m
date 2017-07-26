%--------------------------------------------------------------------------
% bayesfactor.m: Compute numerical estimate of Bayes factor.
%--------------------------------------------------------------------------
%
% DESCRIPTION:
%    Computes numerical estimate of
%
%           Pr(data | H1)
%      BF = ------------- ,
%           Pr(data | H0)
%
%    the probability of the data given the "alternative" hypothesis (H1)
%    over the probability of the data given the "null" hypothesis (H0). This
%    is also known as a Bayes factor (see Kass & Raftery, Journal of the
%    American Statistical Association, 1995). Here we assume that although
%    these probabilities cannot be computed analytically because they
%    involve intractable integrals, we can obtain reasonable estimates of
%    these probabilities with a simple numerical approximation over some
%    latent variable, assuming the prior over the latent variance is
%    uniform. The inputs are the log-probabilities
%
%      Pr(data, Z0 | H0) = Pr(data | Z0, H0) x Pr(Z0 | H0),
%      Pr(data, Z1 | H1) = Pr(data | Z1, H1) x Pr(Z1 | H1),
%
%    where Pr(Z0 | H0) and Pr(Z1 | H1) are uniform over all Z0 and Z1.
%
%    Alternatively, this function can be viewed as computing an importance
%    sampling estimate of the Bayes factor; see, for example, R. M. Neal,
%    "Annealed importance sampling", Statistics and Computing, 2001. This
%    formulation described above is a special case of importance sampling
%    when the settings of the latent variable Z0 and A1 are drawn from the
%    same (uniform) distribution as the prior, Pr(Z0 | H0) and Pr(Z1 | H1),
%    respectively.
%
% USAGE:
%    BF = bayesfactor(logw0, logw1)
%
% INPUT ARGUMENTS:
% logw0  log-probabilities or log-importance weights under H0.
% logw1  log-probabilities or log-importance weights under H1.
%
% OUTPUT ARGUMENTS: The estimated Bayes factor.
%
% LICENSE: GPL v3
%
% DATE: February 3, 2016
%
% AUTHORS:
%    Algorithm was designed by Peter Carbonetto and Matthew Stephens.
%    R, MATLAB and C code was written by Peter Carbonetto.
%    Depts. of Statistics and Human Genetics, University of Chicago,
%    Chicago, IL, USA, and AncestryDNA, San Francisco, CA, USA
%
% REFERENCES:
%    P. Carbonetto, M. Stephens (2012). Scalable variational inference
%    for Bayesian variable selection in regression, and its accuracy in 
%    genetic association studies. Bayesian Analysis 7: 73-108.
%
% SEE ALSO:
%    varbvs, normalizelogweights
%
function BF = bayesfactor (logw0, logw1)

  % Part of the varbvs package, https://github.com/pcarbo/varbvs
  %
  % Copyright (C) 2012-2017, Peter Carbonetto
  %
  % This program is free software: you can redistribute it under the
  % terms of the GNU General Public License; either version 3 of the
  % License, or (at your option) any later version.
  %
  % This program is distributed in the hope that it will be useful, but
  % WITHOUT ANY WARRANY; without even the implied warranty of
  % MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
  % General Public License for more details.
  %
    
  % Compute the marginal log-likelihood under the null hypothesis.
  c     = max(logw0(:));
  logz0 = c + log(mean(exp(logw0(:) - c)));

  % Compute the marginal log-likelihood under the alternative hypothesis.
  c     = max(logw1(:));
  logz1 = c + log(mean(exp(logw1(:) - c)));

  % Compute the numerical estimate of the Bayes factor.
  BF = exp(logz1 - logz0);
