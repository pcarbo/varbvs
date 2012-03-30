\name{updatestats}
\alias{updatestats}
\title{Statistics for updating variational parameters \code{eta}.}

\description{\code{updatestats} calculates useful quantities for
updating the variational approximation to the logistic regression
factors. This function should be called whenever the free parameters
\code{eta} are modified.}

\usage{updatestats(X, y, eta)}

\arguments{

  \item{X}{Matrix of observations about the variables (or features). It
           has \code{n} rows and \code{p} columns, where \code{n} is the
           number of samples, and \code{p} is the number of variables.}

  \item{y}{Vector of observations about the binary outcome. It is a
           vector of length \code{n}.}
	 
  \item{eta}{Free parameters specifying the variational lower bound on
             the logistic regression factors. It is a vector of length
             \code{n}.}}

\value{Returns a list containing five components: \code{xy = t(X) \%*\%
yhat}, \code{xu = t(X) \%*\% u}, \code{d = diag(t(X) \%*\% Uhat \%*\%
X)}, and \code{u = sigmoid(eta - 0.5)/eta}, the slope of the conjugate
to the log-sigmoid function at X, times 2. For a definition of vectors
\code{yhat} and matrix \code{Uhat}, see the Bayesian Analysis paper.}

\author{Peter Carbonetto}

