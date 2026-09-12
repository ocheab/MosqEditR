test_that("tau behaves on simple vectors", {
  expect_equal(mosqedit_tau(c(10, 0, 0)), 1)
  expect_equal(mosqedit_tau(c(1, 1, 1)), 0)
  expect_true(is.na(mosqedit_tau(1)))
})

test_that("percentiles retain NA and rank upward", {
  z <- mosqedit_percentile(c(1, 2, 3, NA))
  expect_equal(z[1:3], c(0, 0.5, 1))
  expect_true(is.na(z[4]))
})

test_that("integrated score uses primary weights", {
  z <- mosqedit_integrated_score(1, 1, 1)
  expect_equal(z$raw_score, 1)
  expect_equal(z$adjusted_score, 1)
  expect_equal(z$domain_coverage, 1)
})

test_that("missing domains are renormalized, not zero-imputed", {
  z <- mosqedit_integrated_score(1, NA, 1)
  expect_equal(z$raw_score, 1)
  expect_equal(z$domain_coverage, 2/3)
  expect_lt(z$adjusted_score, 1)
})

test_that("ranking is deterministic", {
  expect_equal(mosqedit_rank(c(0.8, 0.8, 0.5), c("B", "A", "C")), c(2L, 1L, 3L))
})

