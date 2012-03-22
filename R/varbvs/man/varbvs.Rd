\name{varbvs}
\alias{varbvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varbvs(X, y, sigma, sa, logodds, alpha0 = NULL, mu0 = NULL, verbose = TRUE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
%%     ~~Describe \code{X} here~~
}
  \item{y}{
%%     ~~Describe \code{y} here~~
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
  \item{alpha0}{
%%     ~~Describe \code{alpha0} here~~
}
  \item{mu0}{
%%     ~~Describe \code{mu0} here~~
}
  \item{verbose}{
%%     ~~Describe \code{verbose} here~~
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
function (X, y, sigma, sa, logodds, alpha0 = NULL, mu0 = NULL, 
    verbose = TRUE) 
{
    tolerance <- 1e-04
    if (!is.matrix(X)) 
        stop("Input argument X must be a matrix")
    if (!is.double(X)) 
        X <- double(X)
    n <- nrow(X)
    p <- ncol(X)
    if (length(y) != n) 
        stop("Data X and y do not match")
    if (!is.scalar(sigma) || !is.scalar(sa)) 
        stop("Input arguments sigma and sa must be scalars")
    if (is.scalar(logodds)) 
        logodds <- rep(logodds, p)
    if (length(logodds) != p) 
        stop("Input logodds must be a scalar or a vector of length p")
    if (is.null(alpha0)) {
        alpha <- runif(p)
        alpha <- alpha/sum(alpha)
    }
    else alpha <- alpha0
    if (is.null(mu0)) 
        mu <- rnorm(p)
    else mu <- mu0
    if (length(alpha) != p || length(mu) != p) 
        stop("alpha0 and mu0 must be vectors of length p")
    xy <- c(y \%*\% X)
    d <- diagsq(X)
    Xr <- c(X \%*\% (alpha * mu))
    s <- sa * sigma/(sa * d + 1)
    lnZ <- -Inf
    iter <- 0
    if (verbose) {
        cat("       variational    max. incl max.\n")
        cat("iter   lower bound  change vars E[b]\n")
    }
    while (TRUE) {
        iter <- iter + 1
        alpha0 <- alpha
        mu0 <- mu
        lnZ0 <- lnZ
        params0 <- c(alpha, alpha * mu)
        if (is.odd(iter)) 
            S <- seq(1, p)
        else S <- seq(p, 1, -1)
        result <- varbvsupdate(X, sigma, sa, logodds, xy, d, 
            alpha, mu, Xr, S)
        alpha <- result$alpha
        mu <- result$mu
        Xr <- result$Xr
        lnZ <- intlinear(Xr, d, y, sigma, alpha, mu, s) + intgamma(logodds, 
            alpha) + intklbeta(alpha, mu, s, sigma * sa)
        params <- c(alpha, alpha * mu)
        S <- which(abs(params) > 1e-06)
        err <- relerr(params[S], params0[S])
        if (verbose) 
            cat(sprintf("\%4d \%+13.6e \%0.1e \%4d \%0.2f\n", iter, 
                lnZ, max(err), round(sum(alpha)), max(abs(alpha * 
                  mu))))
        if (lnZ < lnZ0) {
            alpha <- alpha0
            mu <- mu0
            lnZ <- lnZ0
            break
        }
        else if (max(err) < tolerance) 
            break
    }
    return(list(alpha = alpha, mu = mu, s = s, lnZ = lnZ))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
