# Lock tests for the refactored derivative-based threshold methods:
# d2_reg_spline_maxima() (analytic B-spline 2nd derivative, fixed knots) and
# d1_crossing() (Wisen & Wohlfart derivative crossing, front-trimmed, degree 6).
# Golden thresholds guard against regressions in the fits, the analytic
# derivatives, and the shared fit_poly_reg() degree logic.

processed_cpet <- function() {
  cpet <- utils::read.csv(
    system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
  )
  clean <- ventilatory_outliers(
    cpet,
    outlier_cols = "vo2", time = "time", sd_lim = 4, width = 5,
    mos = "mean", align = "center", use_global_sd = TRUE,
    global_sd_mos = "median", exclude_test_val = TRUE,
    remove_outliers = TRUE, max_passes = 1, plot_outliers = FALSE
  )
  avg <- avg_exercise_test(
    clean,
    method = "breath", calc_type = "rolling", time_col = "time", window = 15
  )
  avg$ve_vco2 <- avg$ve / avg$vco2
  avg
}

thr_time <- function(res) {
  t <- res$breakpoint_data$time
  unname(t[!is.na(t)][1])
}

test_that("d2_reg_spline_maxima default (2 knots) is stable at VT2", {
  avg <- processed_cpet()
  res <- d2_reg_spline_maxima(
    avg,
    .x = "time", .y = "ve_vco2", bp = "vt2",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
  )
  expect_true(res$breakpoint_data$determinant_bp)
  expect_equal(thr_time(res), 482.6, tolerance = 1e-2)
})

test_that("d2_reg_spline_maxima knot-search paths run and stay determinant", {
  avg <- processed_cpet()
  for (crit in c("aicc", "lrt")) {
    res <- d2_reg_spline_maxima(
      avg,
      .x = "time", .y = "ve_vco2", bp = "vt2",
      search_knots = TRUE, criterion = crit,
      vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
    )
    expect_true(res$breakpoint_data$determinant_bp)
  }
})

test_that("d1_crossing default (degree 6, front-trimmed) is stable at VT1", {
  avg <- processed_cpet()
  # the crossing is located in time, so .x must be "time"
  res <- d1_crossing(
    avg,
    .x = "time", .y = "vo2", bp = "vt1",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
  )
  expect_equal(thr_time(res), 207.0, tolerance = 1)
})

test_that("d1_crossing degree-search paths run", {
  avg <- processed_cpet()
  for (crit in c("aicc", "lrt")) {
    res <- d1_crossing(
      avg,
      .x = "time", .y = "vo2", bp = "vt1",
      search_degree = TRUE, criterion = crit,
      vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
    )
    expect_false(is.na(thr_time(res)))
  }
})

test_that("shared fit_poly_reg enforces the degree-5 floor", {
  avg <- processed_cpet()
  expect_error(
    d2_poly_reg_maxima(
      avg,
      .x = "time", .y = "ve_vco2", bp = "vt2", degree = 3,
      vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
    ),
    "at least 5"
  )
  expect_error(
    d1_crossing(
      avg,
      .x = "time", .y = "vo2", bp = "vt1", degree = 4,
      vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
    ),
    "at least 5"
  )
})

test_that("invalid criterion is rejected", {
  avg <- processed_cpet()
  expect_error(
    d1_crossing(
      avg,
      .x = "time", .y = "vo2", bp = "vt1", criterion = "bogus",
      vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
    ),
    "criterion"
  )
})

test_that("d2_inflection default (degree 5, front-trimmed) is stable at VT1 and VT2", {
  avg <- processed_cpet()
  vt1 <- d2_inflection(
    avg,
    .x = "vo2", .y = "ve", bp = "vt1",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
  )
  expect_true(vt1$breakpoint_data$determinant_bp)
  expect_equal(unname(vt1$breakpoint_data$vo2), 3826.1, tolerance = 1e-1)

  vt2 <- d2_inflection(
    avg,
    .x = "vco2", .y = "ve", bp = "vt2",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
  )
  expect_true(vt2$breakpoint_data$determinant_bp)
  expect_equal(unname(vt2$breakpoint_data$vco2), 4051.8, tolerance = 1e-1)
})

test_that("d2_inflection honors the ordering argument (was previously ignored)", {
  avg <- processed_cpet()
  by_x <- d2_inflection(
    avg,
    .x = "vo2", .y = "ve", bp = "vt1", ordering = "by_x",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
  )
  by_time <- d2_inflection(
    avg,
    .x = "vo2", .y = "ve", bp = "vt1", ordering = "time",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
  )
  expect_true(by_x$breakpoint_data$determinant_bp)
  expect_true(by_time$breakpoint_data$determinant_bp)
})

test_that("d2_inflection rejects degree < 5 and invalid bp", {
  avg <- processed_cpet()
  expect_error(
    d2_inflection(
      avg,
      .x = "vo2", .y = "ve", bp = "vt1", degree = 3,
      vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
    ),
    "at least 5"
  )
  expect_error(
    d2_inflection(
      avg,
      .x = "vo2", .y = "ve", bp = "vt9",
      vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", plots = FALSE
    ),
    "vt1"
  )
})
