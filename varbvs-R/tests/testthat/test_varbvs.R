context("varbvs")

test_that("model fitting works for simulated data with a continuous outcome",{

  # Run the R script that demonstrates varbvs for mapping of a
  # quantitative trait in a simulated data set.
  source("demo.qtl.R")

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
  expect_equal(summary(fit)$sa$x0,0.158,tolerance = 0.001)

  # Evaluate fitted model.
  expect_equal(cor(y,y.fit)^2,0.650,tolerance = 0.001)
})

test_that(paste("model fitting works for simulated data with a binary",
                "outcome, and no covariates"),{

  # Run the R script that demonstrates mapping of a binary trait in a
  # simulated data set, with no covariates included.
  covariates <- NULL
  source("demo.cc.R",local = TRUE)

  # Check the number of included variables at different probability
  # thresholds.
  expect_equal(summary(fit)$num.included,
               as.table(unlist(list(">0.10" = 6,">0.25" = 6,">0.50" = 5,
                                    ">0.75" = 5,">0.90" = 5,">0.95" = 4))))
  
  # Check the posterior mean of the hyperparameters.
  expect_equal(summary(fit)$logodds$x0,-2.06,tolerance = 0.01)
  expect_equal(summary(fit)$sa$x0,0.209,tolerance = 0.001)
  expect_true(is.na(summary(fit)$sigma$x0))

  # Evaluate fitted model.
  expect_equal(table(factor(y),factor(y.fit)),
               as.table(rbind(c(879,67),
                              c(325,129))),
               check.attributes = FALSE)
})

test_that(paste("model fitting works for simulated data with a binary",
                "outcome, and 2 covariates"),{

  # Run the R script that demonstrates mapping of a binary trait in a
  # simulated data set, with two covariates included in the model.
  covariates <- c("age","weight")
  source("demo.cc.R",local = TRUE)

  # Check the number of included variables at different probability
  # thresholds.
  expect_equal(summary(fit)$num.included,
               as.table(unlist(list(">0.10" = 7,">0.25" = 5,">0.50" = 5,
                                    ">0.75" = 5,">0.90" = 4,">0.95" = 4))))
  
  # Check the posterior mean of the hyperparameters.
  expect_equal(summary(fit)$logodds$x0,-2.07,tolerance = 0.01)
  expect_equal(summary(fit)$sa$x0,0.119,tolerance = 0.001)
  expect_true(is.na(summary(fit)$sigma$x0))

  # Evaluate fitted model.
  expect_equal(table(factor(y),factor(y.fit)),
               as.table(rbind(c(792,124),
                              c(244,240))),
               check.attributes = FALSE)
})

