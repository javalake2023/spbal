# test-spbal-SRS-1.R

testthat::test_that("1. Test SRS.", {
  set.seed(511)
  samp1 <- SRS(total_rows = 10, sample_size = 3)
  expect_equal(samp1, c(5, 1, 10))
})

testthat::test_that("2. Test SRS.", {
  set.seed(511)
  samp2 <- SRS(seed = 69, total_rows = 10, sample_size = 3)
  expect_equal(samp2, c(1, 2, 8))
})

testthat::test_that("3. Test SRS.", {
  set.seed(511)
  samp3 <- SRS(seed = 136, total_rows = 100, sample_size = 4)
  expect_equal(samp3, c(60, 9, 85, 34))
})
