\name{leukemia}
\docType{data}
\alias{leukemia}
\title{Expression levels recorded in leukemia patients.}
\description{
Expression levels recorded for 3,571 genes in 72 patients with leukemia
(Golub et al, 1999). The binary outcome encodes the disease subtype:
acute lymphobastic leukemia (ALL) or acute myeloid leukemia (AML).
}
\usage{data(leukemia)}
\format{
Data are represented as a 72 x 3,571 matrix \code{x} of gene expression
values, and a vector \code{y} of 72 binary disease outcomes.
}
\source{These are the preprocessed data of Dettling (2004) retrieved
from the supplementary materials accompanying Friedman et al (2010).}
\references{
  M. Dettling (2004). BagBoosting for tumor classification with gene
  expression data. \emph{Bioinformatics} \bold{20}, 3583--3593.

  J. Friedman, T. Hastie and R. Tibshirani (2010). Regularization paths
  for generalized linear models via coordinate descent. \emph{Journal of 
  Statistical Software} \bold{33}, 1--22.

  T. R. Golub, et al. (1999). Molecular classification of cancer: class
  discovery and class prediction by gene expression monitoring.
  \emph{Science} \bold{286}, 531--537.
}
\examples{
  # See demo.leukemia.R vignette.
}
\keyword{datasets}
