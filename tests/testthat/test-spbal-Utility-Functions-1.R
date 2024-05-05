# Validate spbal functions, features and parameter validation.

testthat::test_that("1. Verify internal function functions correctly.", {
  sf_object <- sf::st_read(system.file("shape/nc.shp", package="sf"))
  expect_error(spbal::getPanel(sf_object, 1),
               "spbal(getPanel) Simple file object does not contain a feature named panel_id.", fixed=TRUE)
})

# Validate and compare Halton point generators.

#testthat::test_that("2. Compare cppRSHalton() and cppRSHalton_br() same n, seeds.", {
#  chp <- cppRSHalton(n = 1000, seeds = c(123, 456))
#  chpbr <- cppRSHalton_br(n = 1000, seeds = c(123, 456))
#  expect_equal(chp[,2:3], chpbr$pts)
#})

testthat::test_that("3. Compare cppBASpts() and cppRSHalton_br() same n, bases.", {
  chp <- cppBASpts(n = 1000, bases = c(2, 3))
  chpbr <- cppRSHalton_br(n = 1000, bases = c(2, 3), seeds = chp$seeds)
  expect_equal(chp$pts, chpbr$pts)
})

