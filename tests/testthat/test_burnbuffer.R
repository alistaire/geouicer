context("Testing utility functions")
library(geouicer)

test_that("Buffer increments sum to max parameter", {
  expect_equal(sum(buffersteps(1, 5000, 6)), 5000)
})

