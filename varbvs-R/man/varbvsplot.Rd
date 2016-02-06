\name{varbvsplot}
\alias{varbvsplot}
\title{Summarize variable selection results in a single plot.}
\description{
Generate a single plot that summarizes the results of fitting the
Bayesian variable selection model to the data. When the variables are
genetic markers and the groups are chromosomes, the figure resembles a
"Manhattan plot" typically used to summarize the results of a
genome-wide association study or quantitative trait locus (QTL) mapping
study.
}
\usage{
  varbvsplot(fit, score = NULL, vars = NULL, groups = NULL,
             gap = 0, col = "midnightblue",col.vars = "magenta",
             xlab = "", ylab = "",...)
}
\arguments{
  \item{fit}{Output of function \code{varbvs}.}

  \item{cred.int}{Size of credible interval (number between 0 and 1).}

  \item{n}{Show detailed statistics for top n > 0 variables, ranked
  according to their posterior inclusion probabilities.}

  \item{nr}{Number of Monte Carlo samples to draw to estimate credible
    intervals for coefficients of selected variables.}
}
\details{
  Details go here.
}  
\value{
  Input \code{fit} is silently returned.
}
\references{P. Carbonetto and M. Stephens (2012). Scalable variational
inference for Bayesian variable selection in regression, and its
accuracy in genetic association studies. Bayesian Analysis 7(1), pages
73-108.}
\author{Peter Carbonetto <peter.carbonetto@gmail.com>}
\seealso{\code{varbvs}, \code{varbvsprint}, \code{varbvsindep}
}
