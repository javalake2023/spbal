# Validate HIP functions.
# test-spbal-HIP-1.R

testthat::test_that("2. Validate messaging when no population parameter is specified.", {
  # sample size.
  n <- 200
  # number of iterations.
  its <- 7
  # call HIP without population parameter.
  expect_error(spbal::HIP(n = n,
                          iterations =  its), "spbal(HIP) The population parameter must be used. Please specify a population.", fixed=TRUE)
})

