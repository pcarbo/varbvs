\name{varbvsupdate}
\alias{varbvsupdate}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varbvsupdate(X, sigma, sa, logodds, xy, d, alpha0, mu0, Xr0, S)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{sa}{
%%     ~~Describe \code{sa} here~~
}
  \item{logodds}{
%%     ~~Describe \code{logodds} here~~
}
  \item{xy}{
%%     ~~Describe \code{xy} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{alpha0}{
%%     ~~Describe \code{alpha0} here~~
}
  \item{mu0}{
%%     ~~Describe \code{mu0} here~~
}
  \item{Xr0}{
%%     ~~Describe \code{Xr0} here~~
}
  \item{S}{
%%     ~~Describe \code{S} here~~
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
function (X, sigma, sa, logodds, xy, d, alpha0, mu0, Xr0, S) 
{
    if (!is.double(X) || !is.matrix(X)) 
        stop("Input argument X must be a double-precision matrix")
    n <- nrow(X)
    p <- ncol(X)
    m <- length(S)
    if (!is.scalar(sigma) || !is.scalar(sa)) 
        stop("Input arguments sigma and sa must be scalars")
    if (is.scalar(logodds)) 
        logodds <- rep(logodds, p)
    if (length(logodds) != p) 
        stop("Input logodds must be a scalar or a vector of length p")
    if (length(xy) != p || length(d) != p) 
        stop("Inputs xy and d must be vectors of length p")
    if (length(alpha0) != p || length(mu0) != p) 
        stop("Inputs alpha0 and mu0 must be vectors of length p")
    if (length(Xr0) != n) 
        stop("Input Xr0 must be a vector of length n")
    if (sum(S < 1 | S > p) > 0) 
        stop("Input S contains invalid variable indices")
    result <- .C("varbvsupdateR", n = as.integer(n), m = as.integer(m), 
        X = X, sigma = as.double(sigma), sa = as.double(sa), 
        logodds = as.double(logodds), xy = as.double(xy), d = as.double(d), 
        alpha = as.double(alpha0), mu = as.double(mu0), Xr = as.double(Xr0), 
        S = as.integer(S - 1), DUP = FALSE)
    return(list(alpha = result$alpha, mu = result$mu, Xr = result$Xr))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
