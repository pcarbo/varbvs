context("varbvs vignettes")
test_that("leukemia vignette produces correct results",{
  library(knitr)
    
  # Locate the vignette source file, and check that it exists.  
  rmd.file <- system.file("doc/leukemia.Rmd",package = "varbvs")
  expect_gt(nchar(rmd.file),0)

  # Evaluate the code in the leukemia vignette.
  out.file <- tempfile(pattern = "leukemia",fileext = ".R")
  source(purl(rmd.file,output = out.file))

  # Test 1: The number of nonzero coefficients in the fitted glmnet
  # model should be greater than the number of variables included in
  # the varbvs model.
  expect_gt(sum(coef(fit.glmnet,s = lambda.opt) != 0),sum(fit.varbvs$pip))

  # Test 2: The glmnet and varbvs prediction error should be roughly
  # the same (5% tolerance).
  expect_equal(mean(y != y.varbvs),mean(y != y.glmnet),tolerance = 0.05)
})
