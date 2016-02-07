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
