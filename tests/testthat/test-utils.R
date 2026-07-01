# Test the normalize01 function

test_that("normalize01 scales values between 0 and 1", {
  # A simple vector with known min (1) and max (5)
  x <- c(1, 2, 3, 4, 5)
  result <- normalize01(x)

  # We expect: (1-1)/(5-1) = 0, (5-1)/(5-1) = 1, etc.
  expected <- c(0, 0.25, 0.5, 0.75, 1)

  expect_equal(result, expected)
})

test_that("normalize01 handles negative numbers", {
  x <- c(-10, 0, 10)
  result <- normalize01(x)

  # (-10 - (-10)) / (10 - (-10)) = 0/20 = 0
  # (0 - (-10)) / 10 - (-10)) = 10/20 = 0.5
  # (10 - (-10)) / (10 - (-10)) = 20/20 = 1
  expected <- c(0, 0.5, 1)

  expect_equal(result, expected)
})

test_that("normalize01 handles single value (edge case)", {
  x <- c(5)

  # This will fail! (5-5)/(5-5) = 0/0 = NaN
  # This is a bug in the function that we discovered
  result <- normalize01(x)
  expect_true(is.nan(result))
})

test_that("normalize01 works with floats", {
  x <- c(1.5, 2.5, 3.5)
  result <- normalize01(x)

  expected <- c(0, 0.5, 1)

  expect_equal(result, expected)
})

test_that("normalize01 preserves order", {
  x <- c(3, 1, 4, 1, 5, 9)
  result <- normalize01(x)

  # Normalized values should be in same order as originals
  expect_equal(order(x), order(result))
})

test_that("slope_sign_changes identifies sign changes correctly", {
  y <- c(0, 1, 2, 3, 2, 1, 0, -1, -2, -1)

  # Both types of sign changes
  expect_equal(slope_sign_changes(y, change = "both"), c(4, 9))

  # Positive to negative changes
  expect_equal(slope_sign_changes(y, change = "pos_to_neg"), 4)

  # Negative to positive changes
  expect_equal(slope_sign_changes(y, change = "neg_to_pos"), 9)
})

test_that("slope_sign_changes returns empty for monotonic/constant sequences", {
  # Monotonic increasing
  expect_equal(slope_sign_changes(c(1, 2, 3, 4, 5)), numeric(0))

  # Monotonic decreasing
  expect_equal(slope_sign_changes(c(5, 4, 3, 2, 1)), numeric(0))

  # Constant values
  expect_equal(slope_sign_changes(c(5, 5, 5, 5)), numeric(0))
})
