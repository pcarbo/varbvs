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
varbvsplot (fit, score = NULL, groups = NULL, gap = 0,
            col = "midnightblue", vars = NULL, var.labels = NULL,
            var.col = "magenta", pch = 20, xlab = "", ylab = "",
            ltext.args = "col=\"black\",pos=4,cex=0.5", ...)
}
\arguments{
  \item{fit}{Output of function \code{varbvs}.}

  \item{score}{Value to plot on vertical axis. This must be NULL, or a
  numeric vector with one entry for each variable. If score = NULL, the
  posterior probability for each variable is plotted in the vertical
  axis, in which this probability is averaged over hyperparameter
  settings, treating \code{fit$logw} as (unnormalized) log-marginal
  probabilities.}

  \item{groups}{Group the variables in the plot according to this
  argument. This must be equal to NULL, or a vector with one entry for
  each variable. If groups = NULL, all variables are assigned to a
  single group.}

  \item{col}{Argument passed to \code{xyplot} specifying color of points.}

  \item{vars}{Indices of variables to highlight and label.}

  \item{var.labels}{Labels that accompany highlighted variables.}

  \item{var.col}{Color of points for highlighted variables.}

  \item{pch}{Argument passed to \code{xyplot} specifying symbol type.}

  \item{xlab}{Argument passed to \code{xyplot} specifying label for
  horizontal axis.}

  \item{ylab}{Argument passed to \code{xyplot} specifying label for
  vertical axis.}

  \item{ltext.args}{Additional arguments passed to \code{ltext} to draw
  the labels in the plot, specified as a string.}

  \item{...}{Additional arguments passed to \code{xyplot}.}
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
