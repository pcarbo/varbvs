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
})

test_that("model fitting works for simulated data with a binary outcome",{

  # Run the R script that demonstrates mapping of a binary trait in a
  # simulated data set.
  source("demo.cc.R")

  # Check the posterior mean of the hyperparameters.
  expect_equal(summary(fit)$logodds$x0,-2.54,tolerance = 0.01)
  expect_equal(summary(fit)$sa$x0,0.604,tolerance = 0.001)
  expect_true(is.na(summary(fit)$sigma$x0))
}
