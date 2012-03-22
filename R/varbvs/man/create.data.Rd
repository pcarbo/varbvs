\name{create.data}
\alias{create.data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{create.data}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
create.data(snps, sigma, n)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{snps}{
%%     ~~Describe \code{snps} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{n}{
%%     ~~Describe \code{n} here~~
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (snps, sigma, n) 
{
    maf <- snps$maf
    beta <- snps$beta
    p <- length(maf)
    X <- (rnorm(n * p) < rep(maf, n)) + (rnorm(n * p) < rep(maf, 
        n))
    X <- matrix(X, n, p)
    X <- center.columns(X)
    y <- X \%*\% beta + sqrt(sigma) * rnorm(n)
    y <- c(y)
    y <- y - mean(y)
    return(list(X = X, y = y))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
