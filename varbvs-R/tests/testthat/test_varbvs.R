# Part of the varbvs package, https://github.com/pcarbo/varbvs
#
# Copyright (C) 2012-2018, Peter Carbonetto
#
# This program is free software: you can redistribute it under the
# terms of the GNU General Public License; either version 3 of the
# License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANY; without even the implied warranty of
# MERCHANTABILITY OR FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# General Public License for more details.
#
context("varbvs")

test_that("model fitting works for simulated data with a continuous outcome",{

  # Run the varbvs demo for mapping a quantitative trait in a
  # simulated data set.
  demo("varbvs.qtl",package = "varbvs",ask = FALSE)

  # Check the number of included variables at different probability
  # thresholds.
  expect_equal(summary(fit)$num.included,
               as.table(unlist(list(">0.10" = 18,">0.25" = 12,">0.50" = 10,
                                    ">0.75" = 10,">0.90" = 10,">0.95" = 10))))

  # Check the hyperparameter setting with the greatest likelihood value.
  expect_equal(which.max(summary(fit)$logw),10)
  
  # Check the posterior mean of the hyperparameters.
  expect_equal(summary(fit)$logodds$x0,-2.08,tolerance = 0.01)
  expect_equal(summary(fit)$sigma$x0,4.13,tolerance = 0.01)
  expect_equal(summary(fit)$sa$x0,0.158,tolerance = 0.01)

  # Evaluate fitted model.
  expect_equal(cor(y,y.fit)^2,0.650,tolerance = 0.01)
})

test_that(paste("model fitting works for simulated data with a continuous",
                "outcome and one hyperparameter setting"),{

  # The prior inclusion probability is approximately 1/1000.
  logodds <- (-3)
                    
  # Run the varbvs demo for mapping a quantitative trait in a
  # simulated data set.
  source(system.file("demo","varbvs.qtl.R",package = "varbvs"),local = TRUE)

  # Check the posterior mean of the hyperparameters.
  expect_equal(summary(fit)$sigma$x0,4.26,tolerance = 0.01)
  expect_equal(summary(fit)$sa$x0,0.23,tolerance = 0.01)

  # Evaluate fitted model.
  expect_equal(cor(y,y.fit)^2,0.626,tolerance = 0.01)
})

test_that(paste("varbvs model fitting yields accurate estimate of PVE"),{

  # Do not include any additional covariates (Z).
  covariates <- NULL
                    
  # Run the varbvs demo for mapping a quantitative trait in a
  # simulated data set.
  source(system.file("demo","varbvs.qtl.R",package = "varbvs"),local = TRUE)

  # Check that the PVE value used to simulate the data is within the
  # estimated credible interval.
  expect_lte(summary(fit)$model.pve$a,r)
  expect_gte(summary(fit)$model.pve$b,r)
})

test_that(paste("increasing the weights in varbvs results in a",
          "smaller value for sigma"),{

  # Run the varbvs demo for mapping a quantitative trait in a
  # simulated data set.
  demo("varbvs.qtl",package = "varbvs",ask = FALSE)

  # Half the variance of all the observations. This should result in
  # twice as large estimates of sigma.
  fit2 <- varbvs(X,Z,y,"gaussian",logodds = logodds,weights = rep(2,n),n0 = 0)
  expect_equal(range(fit2$sigma/fit$sigma),c(2,2),tolerance = 1e-6)
})

test_that(paste("Bayes factor should support BSLMM in model comparison",
                "demo (varbvs.resid.vcov)"),{

  # Run the model comparison demo.
  demo("varbvs.resid.vcov",package = "varbvs",ask = FALSE)

  # The Bayes factor should provide strong support for the BSLMM model.
  expect_gt(bf,100)
})

test_that(paste("model fitting works for simulated data with a binary",
                "outcome, and no covariates"),{

  # Run the R script that demonstrates mapping of a binary trait in a
  # simulated genetic data set, with no covariates included.
  covariates <- NULL
  source(system.file("demo","varbvs.cc.R",package = "varbvs"),local = TRUE)

  # Check the number of included variables at different probability
  # thresholds.
  expect_equal(summary(fit)$num.included,
               as.table(unlist(list(">0.10" = 6,">0.25" = 6,">0.50" = 5,
                                    ">0.75" = 5,">0.90" = 5,">0.95" = 4))))
  
  # Check the posterior mean of the hyperparameters.
  expect_equal(summary(fit)$logodds$x0,-2.06,tolerance = 0.01)
  expect_equal(summary(fit)$sa$x0,0.209,tolerance = 0.01)
  expect_true(is.na(summary(fit)$sigma$x0))

  # Evaluate fitted model.
  expect_equal(table(factor(y),factor(y.fit)),
               as.table(rbind(c(879,67),
                              c(325,129))),
               check.attributes = FALSE,
               tolerance = 0.02)
})

test_that(paste("model fitting works for simulated data with a binary",
                "outcome, and 2 covariates"),{

  # Run the R script that demonstrates mapping of a binary trait in a
  # simulated genetic data set, with two covariates included.
  covariates <- c("age","weight")
  source(system.file("demo","varbvs.cc.R",package = "varbvs"),local = TRUE)

  # Check the number of included variables at different probability
  # thresholds.
  expect_equal(summary(fit)$num.included,
               as.table(unlist(list(">0.10" = 7,">0.25" = 5,">0.50" = 5,
                                    ">0.75" = 5,">0.90" = 4,">0.95" = 4))))
  
  # Check the posterior mean of the hyperparameters.
  expect_equal(summary(fit)$logodds$x0,-2.07,tolerance = 0.01)
  expect_equal(summary(fit)$sa$x0,0.119,tolerance = 0.01)
  expect_true(is.na(summary(fit)$sigma$x0))

  # Evaluate fitted model.
  expect_equal(table(factor(y),factor(y.fit)),
               as.table(rbind(c(792,124),
                              c(244,240))),
               check.attributes = FALSE,
               tolerance = 0.02)
})

test_that(paste("model fitting works for simulated data with a binary",
                "outcome, 2 covariates, and only 1 hyperparameter setting"),{

  # Run the R script that demonstrates mapping of a binary trait in a
  # simulated genetic data set, with two covariates included.
  logodds    <- (-3)
  covariates <- c("age","weight")
  source(system.file("demo","varbvs.cc.R",package = "varbvs"),local = TRUE)
  
  # Check the posterior mean of the hyperparameters.
  expect_equal(summary(fit)$sa$x0,0.203,tolerance = 0.01)
  expect_true(is.na(summary(fit)$sigma$x0))

  # Evaluate fitted model.
  expect_equal(table(factor(y),factor(y.fit)),
               as.table(rbind(c(794,122),
                              c(248,236))),
               check.attributes = FALSE,
               tolerance = 0.02)
})

test_that("model fitting works when crossprod(Z) is near-singular",{

  # Load the data with a near-singular covariate matrix Z. To see that
  # crossprod(Z) is near-singular, run this code:
  #
  #  R <- cor(Z)
  #  v <- eigen(R)$values
  #  print(v)
  #
  load(system.file("datafiles","singular.RData",package = "varbvs"))
  expect_silent(fit <- varbvs(X,Z,y,"gaussian",verbose = FALSE))
})

test_that(paste("model fitting works for linear regression with",
                "mixture-of-normals priors in simulated data"),{
  tol <- 0.05
                    
  # Run the R script that demonstrates varbvsmix on a simulated data
  # set in which all the candidate variables are uncorrelated.
  demo("varbvsmix",package = "varbvs",ask = FALSE)

  # The variational lower bound should always be increasing.
  expect_true(all(diff(fit$logZ) > 0))
  
  # Check the estimates of the model parameters against the
  # settings used to simulate the data.
  expect_equal(w,fit$w,tolerance = tol)
  expect_equal(fit$sigma,se,tolerance = tol*se)
})

test_that("varbvs and varbvsmix produce same estimates when K=2",{
  tol <- 1e-4
  demo("varbvsmix.test",package = "varbvs",ask = FALSE)
  rownames(fit$alpha) <- NULL
  rownames(fit$mu)    <- NULL
  
  # Check that the varbvs and varbvsmix parameter estimates are the
  # same.
  niter <- length(fit$logZ)
  expect_equal(fit$alpha[,2],c(fit2$alpha),tolerance = tol)
  expect_equal(fit$mu[,2],c(fit2$mu),tolerance = tol)
  expect_equal(fit$logZ[niter],fit2$logw,tolerance = tol)
})
