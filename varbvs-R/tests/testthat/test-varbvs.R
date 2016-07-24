context("varbvs")

test_that("model fitting works for continuous outcome in simulated data",{
  source("demo.qtl.R")
  expect_equal(summary(fit)$num.included,
               as.table(unlist(list(">0.10" = 18,">0.25" = 12,">0.50" = 10,
                                    ">0.75" = 10,">0.90" = 10,">0.95" = 10))))
})

