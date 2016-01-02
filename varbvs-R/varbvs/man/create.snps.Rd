\name{create.snps}
\alias{create.snps}
\title{Generate SNP Info}

\description{Generates minor allele frequencies and additive effects for
genetic loci (specifically, these are single nucleotide polymorphisms,
or SNPs for short).}

\usage{create.snps(p, n)}

\arguments{
  \item{p}{number of SNPs.}
  \item{n}{number of causal SNPs.}
}

\details{Additive effects are generated from the standard normal, and
minor allele frequencies are uniform between 0.05 and 0.5.}

\value{
Returns a list containing two components:
\item{maf}{vector with minor allele frequencies of SNPs.}
\item{beta}{vector with additive effects of SNPs.}
}

\seealso{
\code{\link{create.data}}
}

\examples{

## Generate minor allele frequencies and additive effects for 1000
## SNPs; of these SNPs, 10 have nonzero additive effects on the trait.
snps <- create.snps(1000,10)

}

\author{Peter Carbonetto}
