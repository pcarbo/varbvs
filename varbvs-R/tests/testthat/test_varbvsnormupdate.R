context("varbvsnormupdate")
test_that("all versions of varbvsnormupdate produce the same result",{

  # SCRIPT PARAMETERS
  # -----------------
  n       <- 2400  # Number of samples.
  p       <- 5000  # Number of variables (genetic markers).
  na      <- 50    # Number of quantitative trait loci (QTLs).
  se      <- 4     # Variance of residual.
  r       <- 0.5   # Proportion of variance in trait explained by QTLs.
  logodds <- -2    # Prior log-odds of inclusion.

  # Set the random number generator seed.
  set.seed(1)

  # GENERATE DATA SET
  # -----------------
  # Generate the minor allele frequencies so that they are uniform
  # over range [0.05,0.5]. Then simulate genotypes assuming all
  # markers are uncorrelated (i.e., unlinked), according to the
  # specified minor allele frequencies.
  cat("Generating data set.\n")
  maf <- 0.05 + 0.45*runif(p)
  X   <- (runif(n*p) < maf) +
         (runif(n*p) < maf)
  X   <- matrix(as.double(X),n,p,byrow = TRUE)

  # Generate additive effects for the markers so that exactly na of
  # them have a nonzero effect on the trait.
  i       <- sample(p,na)
  beta    <- rep(0,p)
  beta[i] <- rnorm(na)

  # Adjust the QTL effects so that we control for the proportion of variance
  # explained (r). That is, we adjust beta so that r = a/(a+1), where I've
  # defined a = beta'*cov(X)*beta. Here, sb is the variance of the (nonzero)
  # QTL effects.
  sb   <- r/(1-r)/var(c(X %*% beta))
  beta <- sqrt(sb*se) * beta

  # Generate the quantitative trait measurements.
  y <- c(X %*% beta + sqrt(se)*rnorm(n))

  # Center the columns of X and subtracting the mean from y.
  rep.row <- function (x, n)
    matrix(x,n,length(x),byrow = TRUE)
  X <- X - rep.row(colMeans(X),n)
  y <- y - mean(y)

  # RUN CO-ORDINATE ASCENT UPDATES
  # ------------------------------
  cat("Running one round of co-ordinate ascent updates.\n")

  # Randomly set initial estimates of the variational parameters.
  alpha0 <- runif(p)
  alpha0 <- alpha0/sum(alpha0)
  mu0    <- rnorm(p)
  
  # Set up the inputs to varbvsnormupdate.
  xy      <- c(y %*% X)
  d       <- diagsq(X)
  Xr0     <- c(X %*% (alpha0*mu0))
  logodds <- rep(log(10)*logodds,p)

  r <- cbind(system.time(out1 <- varbvsnormupdate(X,se,sb,logodds,xy,d,
                                                  alpha0,mu0,Xr0,1:p,"R")),
             system.time(out2 <- varbvsnormupdate(X,se,sb,logodds,xy,d,
                                                  alpha0,mu0,Xr0,1:p,".Call")),
             system.time(out3 <- varbvsnormupdate(X,se,sb,logodds,xy,d,
                                                  alpha0,mu0,Xr0,1:p,"Rcpp")))
  r        <- as.data.frame(r)
  names(r) <- c("R",".Call","Rcpp")
  print(r["elapsed",])

  # Check the outputs from all three versions of varbvsnormupdate.
  expect_equal(out1$alpha,out2$alpha,tolerance = 1e-8)
  expect_equal(out1$alpha,out3$alpha,tolerance = 1e-8)
  expect_equal(out1$mu,out2$mu,tolerance = 1e-8)
  expect_equal(out1$mu,out3$mu,tolerance = 1e-8)
  expect_equal(out1$Xr,out2$Xr,tolerance = 1e-8)
  expect_equal(out1$Xr,out3$Xr,tolerance = 1e-8)
})
