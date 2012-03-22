\name{intlinear}
\alias{intlinear}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
intlinear(Xr, d, y, sigma, alpha, mu, s)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{Xr}{
%%     ~~Describe \code{Xr} here~~
}
  \item{d}{
%%     ~~Describe \code{d} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
}
  \item{sigma}{
%%     ~~Describe \code{sigma} here~~
}
  \item{alpha}{
%%     ~~Describe \code{alpha} here~~
}
  \item{mu}{
%%     ~~Describe \code{mu} here~~
}
  \item{s}{
%%     ~~Describe \code{s} here~~
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
function (Xr, d, y, sigma, alpha, mu, s) 
{
    n <- length(y)
    f <- -n/2 * log(2 * pi * sigma) - norm2(y - Xr)^2/(2 * sigma) - 
        dot(d, betavar(alpha, mu, s))/(2 * sigma)
    return(f)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
