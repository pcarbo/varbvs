\name{varsimbvs}
\alias{varsimbvs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
varsimbvs(X, y, sigma, sa, log10q, a, b, ca)
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
  \item{log10q}{
%%     ~~Describe \code{log10q} here~~
}
  \item{a}{
%%     ~~Describe \code{a} here~~
}
  \item{b}{
%%     ~~Describe \code{b} here~~
}
  \item{ca}{
%%     ~~Describe \code{ca} here~~
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
function (X, y, sigma, sa, log10q, a, b, ca) 
{
    as <- 0.01
    bs <- 0.01
    if (!is.matrix(X)) 
        stop("Input argument X must be a matrix")
    n <- nrow(X)
    p <- ncol(X)
    ns <- length(sigma)
    if (!is.scalar(a) || !is.scalar(b) || !is.scalar(c)) 
        stop("Input arguments a, b and c must be scalars")
    if (length(sa) != ns || length(log10q) != ns) 
        stop("Inputs sigma, sa and log10q must be the same size")
    sx <- sum(sd(X)^2)
    q1 <- 10^log10q
    lnZ <- rep(0, ns)
    logw <- rep(0, ns)
    alpha <- matrix(0, p, ns)
    mu <- matrix(0, p, ns)
    cat(sprintf("Finding best initialization for \%d combinations ", 
        ns))
    cat("of hyperparameters.\n")
    for (i in 1:ns) {
        cat(sprintf("(\%03d) sigma = \%4.1f, sa = \%0.3f, q = \%0.2e", 
            i, sigma[i], sa[i], q1[i]))
        cat(rep("\b", 1, 44), sep = "")
        alpha0 <- runif(p)
        alpha0 <- alpha0/sum(alpha0)
        mu0 <- rnorm(p)
        result <- varbvs(X, y, sigma[i], sa[i], logit(q1[i]), 
            alpha0, mu0, verbose = FALSE)
        lnZ[i] <- result$lnZ
        alpha[, i] <- result$alpha
        mu[, i] <- result$mu
    }
    cat("\n")
    i <- which.max(lnZ)
    alpha0 <- c(alpha[, i])
    mu0 <- c(mu[, i])
    cat(sprintf("Computing importance weights for \%d combinations ", 
        ns))
    cat("of hyperparameters.\n")
    for (i in 1:ns) {
        cat(sprintf("(\%03d) sigma = \%4.1f, sa = \%0.3f, q = \%0.2e", 
            i, sigma[i], sa[i], q1[i]))
        cat(rep("\b", 1, 44), sep = "")
        result <- varbvs(X, y, sigma[i], sa[i], logit(q1[i]), 
            alpha0, mu0, verbose = FALSE)
        logw[i] <- lnZ[i] + loginvgamma(sigma[i], as/2, bs/2) + 
            logpve(ca * sx, sa[i]) + dbeta(q1[i], a + 1, b, log = TRUE)
    }
    cat("\n")
    w <- normalizelogweights(logw)
    alpha <- c(alpha \%*\% w)
    mu <- c(mu \%*\% w)
    return(list(w = w, alpha = alpha, mu = mu))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
