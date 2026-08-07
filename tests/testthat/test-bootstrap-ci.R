# Tests for the shared bootstrap CI engine (bootstrap_threshold_ci()) and its
# use across all four derivative-based threshold methods, plus the exact
# polynomial-arithmetic helpers (poly_deriv_coefs/poly_eval_coefs) that
# replaced symbolic differentiation + Ryacas in the hot path of
# d1_crossing()/d2_inflection().

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
  avg_exercise_test(
    clean,
    method = "breath", calc_type = "rolling", time_col = "time", window = 15
  )
}

# ---------------------------------------------------------------------------
# poly_deriv_coefs() / poly_eval_coefs(): exact arithmetic vs. symbolic Deriv
# ---------------------------------------------------------------------------

test_that("poly_deriv_coefs matches symbolic Deriv::Deriv", {
  coefs <- c(2, -3, 0.5, 1.2, -0.01, 0.0003) # arbitrary degree-5 poly
  poly_expr <- expr_from_coefs(coefs)

  deriv1_sym <- expr_to_func(Deriv::Deriv(poly_expr, x = "x", nderiv = 1))
  deriv2_sym <- expr_to_func(Deriv::Deriv(poly_expr, x = "x", nderiv = 2))

  deriv1_coefs <- poly_deriv_coefs(coefs)
  deriv2_coefs <- poly_deriv_coefs(deriv1_coefs)

  xe <- seq(-10, 10, length.out = 100)
  expect_equal(poly_eval_coefs(deriv1_coefs, xe), deriv1_sym(xe), tolerance = 1e-8)
  expect_equal(poly_eval_coefs(deriv2_coefs, xe), deriv2_sym(xe), tolerance = 1e-8)
})

test_that("poly_deriv_coefs handles a constant (derivative is 0)", {
  expect_equal(poly_deriv_coefs(5), 0)
})

test_that("poly_subtract_coefs pads the shorter vector", {
  expect_equal(poly_subtract_coefs(c(1, 2, 3), c(1, 1)), c(0, 1, 3))
})

# ---------------------------------------------------------------------------
# bootstrap_threshold_ci(): NA-handling engine
# ---------------------------------------------------------------------------

test_that("bootstrap_threshold_ci computes a sane interval when replicates are clean", {
  set.seed(1)
  d <- data.frame(x = 1:50)
  # every replicate succeeds and returns a value near 10
  always_valid <- function(.data, i) mean(.data$x[i]) # ~25.5, deterministic-ish
  res <- bootstrap_threshold_ci(d, statistic = always_valid, R = 199, conf_level = 0.95)
  expect_false(is.null(res))
  expect_equal(res$n_valid, res$n_total)
  expect_true(res$lower < res$upper)
})

test_that("bootstrap_threshold_ci warns and excludes NAs when some replicates fail", {
  set.seed(1)
  d <- data.frame(x = 1:50)
  # ~20% of replicates fail
  sometimes_na <- function(.data, i) {
    if (stats::runif(1) < 0.2) {
      return(NA_real_)
    }
    mean(.data$x[i])
  }
  expect_warning(
    res <- bootstrap_threshold_ci(d, statistic = sometimes_na, R = 199, conf_level = 0.95),
    "did not produce a valid threshold"
  )
  expect_false(is.null(res))
  expect_true(res$n_valid < res$n_total)
  expect_true(res$lower < res$upper)
})

test_that("bootstrap_threshold_ci returns NULL and warns when almost all replicates fail", {
  set.seed(1)
  d <- data.frame(x = 1:50)
  almost_always_na <- function(.data, i) NA_real_
  expect_warning(
    res <- bootstrap_threshold_ci(d, statistic = almost_always_na, R = 199, conf_level = 0.95),
    "could not be reliably estimated"
  )
  expect_null(res)
})

# ---------------------------------------------------------------------------
# End-to-end: ci = TRUE on all four derivative-based methods
# ---------------------------------------------------------------------------

test_that("all four derivative methods produce ordered, correctly-labeled CIs", {
  avg <- processed_cpet()
  avg$ve_vco2 <- avg$ve / avg$vco2

  check_ci <- function(res, algo_name) {
    d <- res$breakpoint_data
    expect_equal(nrow(d), 3)
    expect_setequal(d$est_ci, c("lower_ci", "estimate", "upper_ci"))
    expect_true(all(d$determinant_bp))
    expect_true(all(d$algorithm == algo_name))
    x_col <- setdiff(names(d), c(
      "bp", "algorithm", "x_var", "y_var", "est_ci", "determinant_bp"
    ))
    x_col
  }

  set.seed(1)
  r1 <- d2_poly_reg_maxima(
    avg,
    .x = "time", .y = "ve_vco2", bp = "vt2",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", ci = TRUE, plots = FALSE
  )
  check_ci(r1, "d2_poly_reg_maxima")
  d1 <- r1$breakpoint_data
  expect_true(d1$time[d1$est_ci == "lower_ci"] < d1$time[d1$est_ci == "estimate"])
  expect_true(d1$time[d1$est_ci == "estimate"] < d1$time[d1$est_ci == "upper_ci"])

  set.seed(1)
  r2 <- d2_reg_spline_maxima(
    avg,
    .x = "time", .y = "ve_vco2", bp = "vt2",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", ci = TRUE, plots = FALSE
  )
  check_ci(r2, "d2_reg_spline_maxima")

  set.seed(1)
  r3 <- d1_crossing(
    avg,
    .x = "time", .y = "vo2", bp = "vt1",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", ci = TRUE, plots = FALSE
  )
  check_ci(r3, "d1_crossing")

  set.seed(1)
  r4 <- d2_inflection(
    avg,
    .x = "vo2", .y = "ve", bp = "vt1",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", ci = TRUE, plots = FALSE
  )
  check_ci(r4, "d2_inflection")
})

test_that("ci = TRUE is skipped cleanly when the point estimate is indeterminate", {
  cpet <- utils::read.csv(
    system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
  )
  cpet_20s <- avg_exercise_test(
    cpet,
    method = "time", calc_type = "bin", time_col = "time", window = 20
  )
  cpet_20s$ve_vco2 <- cpet_20s$ve / cpet_20s$vco2

  res <- d2_poly_reg_maxima(
    cpet_20s,
    .x = "time", .y = "ve_vco2", bp = "vt2",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", ci = TRUE, plots = FALSE
  )
  # only the indeterminate estimate row; no bootstrap CI rows were appended
  expect_equal(nrow(res$breakpoint_data), 1)
  expect_false(res$breakpoint_data$determinant_bp)
})

test_that("d2_poly_reg_maxima CI rows are labeled with the correct algorithm", {
  # regression test for a copy-paste bug where CI rows were mislabeled
  # "d2_reg_spline_maxima"
  avg <- processed_cpet()
  avg$ve_vco2 <- avg$ve / avg$vco2
  set.seed(1)
  res <- d2_poly_reg_maxima(
    avg,
    .x = "time", .y = "ve_vco2", bp = "vt2",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time", ci = TRUE, plots = FALSE
  )
  expect_true(all(res$breakpoint_data$algorithm == "d2_poly_reg_maxima"))
})
