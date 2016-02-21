# Returns the prior variance of the regression coefficients given the
# prior estimate of the proportion of variance explained by the
# variables (h), and the prior log-odds that a variable is included in
# the model (logodds). Note that the base-10 logarithm is used here to
# define the prior log-odds.
pve2sa <- function (X, h, logodds) {
  sx <- sum(var1.cols(X))
  return(h/(1-h)/(sigmoid(log(10)*logodds)*sx))
}
