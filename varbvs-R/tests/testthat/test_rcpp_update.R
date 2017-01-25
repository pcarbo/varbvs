context("Rcpp Updates")

test_that("Rcpp Update Works", {
  ## generate test data ------------------------------------------------------
  gen_dat <- function(seed) {
    set.seed(seed)
    n <- 11
    p <- 21
    X <- matrix(stats::rnorm(n * p), nrow = n)
    beta <- stats::rnorm(p)
    sigma <- 1
    sa <- 1.1
    logodds <- rep(1, p)
    y <- X %*% beta + stats::rnorm(n = n, sd = sqrt(sigma))
    xy <- crossprod(X, y)
    d <- colSums(X ^ 2)
    alpha <- stats::runif(n = p)
    mu <- rep(0, p)
    Xr <- X %*% (mu * alpha)
    i <- sample(1:p)
    return(list(X = X, sigma = sigma, sa = sa, logodds = logodds, xy = xy,
                d = d, alpha = alpha, mu = mu, Xr = Xr, i = i))
  }

  ## cpp version -------------------------------------------------------------
  dat1 <- gen_dat(101)
  varbvsnormupdate_rcpp(X = dat1$X, sigma = dat1$sigma, sa = dat1$sa,
                        logodds = dat1$logodds, xy = dat1$xy,
                        d = dat1$d, alpha = dat1$alpha, mu = dat1$mu,
                        Xr = dat1$Xr, i = dat1$i - 1)

  ## .Call version -----------------------------------------------------------
  dat2 <- gen_dat(101)
  out <- .Call("varbvsnormupdate_Call", X = dat2$X,
               sigma = as.double(dat2$sigma),
               sa = as.double(dat2$sa), logodds = as.double(dat2$logodds),
               xy = as.double(dat2$xy), d = as.double(dat2$d),
               alpha = dat2$alpha, mu = dat2$mu, Xr = dat2$Xr,
               i = as.integer(dat2$i - 1))

  ## R version ---------------------------------------------------------------
  dat3 <- gen_dat(101)
  for (j in dat3$i) {

    # Compute the variational estimate of the posterior variance.
    s <- dat3$sa * dat3$sigma / (dat3$sa * dat3$d[j] + 1)

    # Update the variational estimate of the posterior mean.
    r     <- dat3$alpha[j] * dat3$mu[j]

    dat3$mu[j] <- s / dat3$sigma * (dat3$xy[j] + dat3$d[j] * r -
                                      sum(dat3$X[, j] * dat3$Xr))

    # Update the variational estimate of the posterior inclusion
    # probability.
    dat3$alpha[j] <- sigmoid(dat3$logodds[j] +
                               (log(s / (dat3$sa * dat3$sigma)) +
                                  dat3$mu[j] ^ 2 / s) / 2)

    # Update Xr = X*r.
    dat3$Xr <- dat3$Xr + (dat3$alpha[j] * dat3$mu[j] - r) * dat3$X[,j]
  }

  ## Test --------------------------------------------------------------------
  expect_equal(dat1$mu, dat2$mu, dat3$mu)
  expect_equal(dat1$alpha, dat2$alpha, dat3$alpha)
  expect_equal(dat1$Xr, dat2$Xr, dat3$Xr)

}
)

test_that("interface with varbvsnormupdate is ok", {
  ## generate test data ------------------------------------------------------
  gen_dat <- function(seed) {
    set.seed(seed)
    n <- 11
    p <- 21
    X <- matrix(stats::rnorm(n * p), nrow = n)
    beta <- stats::rnorm(p)
    sigma <- 1
    sa <- 1.1
    logodds <- rep(1, p)
    y <- X %*% beta + stats::rnorm(n = n, sd = sqrt(sigma))
    xy <- crossprod(X, y)
    d <- colSums(X ^ 2)
    alpha <- stats::runif(n = p)
    mu <- rep(0, p)
    Xr <- X %*% (mu * alpha)
    i <- sample(1:p)
    return(list(X = X, sigma = sigma, sa = sa, logodds = logodds, xy = xy,
                d = d, alpha = alpha, mu = mu, Xr = Xr, i = i))
  }

  dat <- gen_dat(2345234)
  out1 <- varbvsnormupdate(X = dat$X, sigma = dat$sigma, sa = dat$sa,
                           logodds = dat$logodds, xy = dat$xy, d = dat$d,
                           alpha0 = dat$alpha, mu0 = dat$mu,
                           Xr0 = dat$Xr, i = dat$i,algorithm.version = "Rcpp")
  out2 <- varbvsnormupdate(X = dat$X, sigma = dat$sigma, sa = dat$sa,
                           logodds = dat$logodds, xy = dat$xy, d = dat$d,
                           alpha0 = dat$alpha, mu0 = dat$mu,
                           Xr0 = dat$Xr, i = dat$i,
                           algorithm.version = ".Call")
  out3 <- varbvsnormupdate(X = dat$X, sigma = dat$sigma, sa = dat$sa,
                           logodds = dat$logodds, xy = dat$xy, d = dat$d,
                           alpha0 = dat$alpha, mu0 = dat$mu,
                           Xr0 = dat$Xr, i = dat$i,
                           algorithm.version = "R")

  expect_equal(out1$alpha, out2$alpha, out3$alpha)
  expect_equal(out1$mu, out2$mu, out3$mu)
  expect_equal(out1$Xr, out2$Xr, out3$Xr)
}
)
