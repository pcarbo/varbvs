\name{subset.varbvs}
\alias{subset.varbvs}
\title{Select hyperparameter settings from varbvs analysis.}
\description{Select a subset of the candidate hyperparameter settings,
and return a new varbvs analysis object with these hyperparameter
settings only.}
\usage{
  \method{subset}{varbvs}(x, subset, \ldots)
}
\arguments{
  \item{x}{Output of function \code{\link{varbvs}}.}

  \item{subset}{Expression indicating hyperparameter settings to select.
    This expression should include one or more of \code{logodds},
    \code{sigma} and \code{sa}.}

  \item{\dots}{Other arguments to generic subset function. These
    extra arguments are not used here.}
}
\value{
  An object with S3 class \code{c("varbvs","list")}.
}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. \emph{Bayesian Analysis}
\bold{7}, 73--108.}
\author{Peter Carbonetto \email{peter.carbonetto@gmail.com}}
\seealso{\code{varbvs}}
\examples{

  # First run one of the examples in help(varbvs), then try running 
  # this.
  #
  #   fit.new <- subset(fit,logodds < (-2))
  # 
}
