# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2017, Peter Carbonetto
#
# This program is free software: you can redistribute it under the
# terms of the GNU General Public License; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
# Compute numerical estimate of Bayes factor.
bayesfactor <- function (logw0, logw1) {

  # Compute the marginal log-likelihood under the null hypothesis.
  c     <- max(logw0)
  logz0 <- c + log(mean(exp(logw0 - c)))

  # Compute the marginal log-likelihood under the alternative hypothesis.
  c     <- max(logw1)
  logz1 <- c + log(mean(exp(logw1 - c)))

  # Compute the numerical estimate of the Bayes factor.
  return(exp(logz1 - logz0))
}
