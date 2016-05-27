\name{cytokine}
\docType{data}
\alias{cytokine}
\title{Cytokine signaling genes SNP annotation.}
\description{This gene set was selected in Carbonetto and Stephens
 (2013) from an interrogation of 3,158 derived from 8 publicly available
 pathway databases.
}
\usage{data(cytokine)}
\format{
\code{cytokine[i] = 1} if SNP i lies within 100 kb of a gene in the "Cytokine
signaling in immune system" gene set, and \code{cytokine[i] = 0} otherwise.
}
\source{Pathway id 75790 from the Reactome database, or pathway id
  366171 from the BioSystems database.}
\references{
  P. Carbonetto and M. Stephens (2013). Integrated enrichment analysis of
  variants and pathways in genome-wide association studies indicates
  central role for IL-2 signaling genes in type 1 diabetes, and cytokine
  signaling genes in Crohn's disease. \emph{PLoS Genetics} \bold{9},
  e1003770.
}
\examples{
  # See demo.cytokine.R vignette.
}
\keyword{datasets}
