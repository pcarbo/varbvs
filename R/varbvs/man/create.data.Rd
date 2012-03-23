\name{create.data}
\alias{create.data}
\title{Simulate SNP data}

\description{Generates samples of the genotypes and quantitative trait
  (i.e. a continuous outcome Y) according to the specified SNP minor
  allele frequencies and additive effects.}

\usage{create.data(maf, beta, sigma, n)}

\arguments{
  \item{maf}{vector giving minor allele frequencies of the SNPs.}
  \item{beta}{vector of the same length giving the additive effects
              of the SNPs.}
  \item{sigma}{positive scalar specifying the variance of residual.}
  \item{n}{number of samples to generate.}
}

\details{Genotypes are generated from a binomial distribution with
  success rates given by the minor allele frequencies. Observations
  about the quantitative trait are generated according to \code{y =
  X*beta + e}, where the residual \code{e} is normal with mean zero and
  covariance \code{sigma} times the identity matrix.}

\value{
Returns a list containing two components:
\item{X}{matrix of genotype data with centered columns so that the mean
         of each column is zero.}
\item{y}{vector of quantitative trait data, centered so that the mean of
         the vector entries is zero.}
}

\author{Peter Carbonetto}

\seealso{\code{\link{create.snps}}}

\examples{

## Randomly generate genotypes and quantitative trait measurements for
## 1000 SNPs and 200 individuals, in which the variance of the residual
## error is 4. Of these SNPs, 10 have nonzero additive effects on the
## trait.
se   <- 4
snps <- create.snps(1000,10)
data <- create.data(snps$maf,snps$beta,se,200)

## Calculate the proportion of variance in the quantitative trait
## explained by the SNPs.
sz <- var(c(data$X \%*\% snps$beta))
sz/(sz + se)

}