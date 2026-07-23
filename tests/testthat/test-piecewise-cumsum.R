# Lock the cumulative-sufficient-statistics search engines
# (piecewise_metrics_cumsum / piecewise_metrics_cumsum_jm) to an independent
# brute-force reference: fitting a fresh pair of lm() models at every split.
# These oracles reimplement the search the slow way, so the fast engines can
# never silently drift from the published algorithms.

# --- brute-force oracles (lm per split) -------------------------------------

# ORR / v-slope flavor: two unconstrained OLS lines sharing the split row.
brute_unconstrained <- function(d, .x, .y, conf_level = 0.95) {
  x <- d[[.x]]
  y <- d[[.y]]
  n <- length(x)
  RSS_two <- MSE_two <- pf_two <- int_point_x <- rep(NA_real_, n)
  pos_change <- pos_slope_after_bp <- rep(NA, n)

  lm_s <- stats::lm(y ~ x)
  RSS_s <- sum(stats::resid(lm_s)^2)

  for (i in 2:(n - 1)) {
    L <- stats::lm(y[1:i] ~ x[1:i])
    R <- stats::lm(y[i:n] ~ x[i:n]) # shared split row
    RSS_two[i] <- sum(stats::resid(L)^2) + sum(stats::resid(R)^2)
    MSE_two[i] <- RSS_two[i] / (n - 4)
    f <- (RSS_s - RSS_two[i]) / (2 * MSE_two[i])
    pf_two[i] <- stats::pf(f, 2, n - 4, lower.tail = FALSE)
    bL <- stats::coef(L)
    bR <- stats::coef(R)
    pos_change[i] <- 100 * (bR[2] - bL[2]) / abs(bL[2]) > 0
    pos_slope_after_bp[i] <- bR[2] > 0
    int_point_x[i] <- (bR[1] - bL[1]) / (bL[2] - bR[2])
  }

  crit_F <- stats::qf(conf_level, 1, n - 4)
  inside_ci <- (RSS_two - min(RSS_two, na.rm = TRUE)) /
    min(MSE_two, na.rm = TRUE) < crit_F
  tibble::tibble(
    p = pf_two, pos_change = pos_change,
    pos_slope_after_bp = pos_slope_after_bp,
    int_point_x = int_point_x, inside_ci = inside_ci, idx = seq_len(n)
  )
}

# JM flavor: left OLS, right regression constrained through the join point,
# right segment does NOT share the split row.
brute_jm <- function(d, .x, .y, conf_level = 0.95) {
  x <- d[[.x]]
  y <- d[[.y]]
  n <- length(x)
  RSS_two <- MSE_two <- pf_two <- int_point_x <- rep(NA_real_, n)
  pos_change <- pos_slope_after_bp <- rep(NA, n)

  lm_s <- stats::lm(y ~ x)
  RSS_s <- sum(stats::resid(lm_s)^2)

  for (i in 2:(n - 2)) { # exclude 1, n-1, n
    L <- stats::lm(y[1:i] ~ x[1:i])
    bL <- stats::coef(L)
    x0 <- x[i]
    y0 <- bL[1] + bL[2] * x0
    xr <- x[(i + 1):n] # strictly greater than x0
    yr <- y[(i + 1):n]
    u <- xr - x0
    v <- yr - y0
    b3 <- sum(u * v) / sum(u^2) # through-origin slope
    RSS_right <- sum((v - b3 * u)^2)
    RSS_two[i] <- sum(stats::resid(L)^2) + RSS_right
    MSE_two[i] <- RSS_two[i] / (n - 4)
    f <- (RSS_s - RSS_two[i]) / (2 * MSE_two[i])
    pf_two[i] <- stats::pf(f, 2, n - 4, lower.tail = FALSE)
    pos_change[i] <- 100 * (b3 - bL[2]) / abs(bL[2]) > 0
    pos_slope_after_bp[i] <- b3 > 0
    int_point_x[i] <- x0 # JM breakpoint is the observed x
  }

  crit_F <- stats::qf(conf_level, 1, n - 4)
  inside_ci <- (RSS_two - min(RSS_two, na.rm = TRUE)) /
    min(MSE_two, na.rm = TRUE) < crit_F
  tibble::tibble(
    p = pf_two, pos_change = pos_change,
    pos_slope_after_bp = pos_slope_after_bp,
    int_point_x = int_point_x, inside_ci = inside_ci, idx = seq_len(n)
  )
}

# small, quick dataset for the oracle comparisons
test_data <- function() {
  cpet <- utils::read.csv(
    system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
  )
  avg_exercise_test(
    cpet,
    method = "time", calc_type = "bin", time_col = "time", window = 30
  )
}

# --- tests ------------------------------------------------------------------

test_that("piecewise_metrics_cumsum matches brute-force OLS (ORR/v-slope)", {
  d <- test_data()
  cs <- piecewise_metrics_cumsum(d, .x = "vco2", .y = "ve")
  bf <- brute_unconstrained(d, .x = "vco2", .y = "ve")

  expect_equal(cs$p, bf$p, tolerance = 1e-6)
  expect_equal(cs$int_point_x, bf$int_point_x, tolerance = 1e-6)
  expect_identical(cs$pos_change, bf$pos_change)
  expect_identical(cs$pos_slope_after_bp, bf$pos_slope_after_bp)
  expect_identical(cs$inside_ci, bf$inside_ci)
})

test_that("piecewise_metrics_cumsum_jm matches brute-force constrained JM", {
  d <- test_data()
  cs <- piecewise_metrics_cumsum_jm(d, .x = "vco2", .y = "ve")
  bf <- brute_jm(d, .x = "vco2", .y = "ve")

  expect_equal(cs$p, bf$p, tolerance = 1e-6)
  expect_equal(cs$int_point_x, bf$int_point_x, tolerance = 1e-6)
  expect_identical(cs$pos_change, bf$pos_change)
  expect_identical(cs$pos_slope_after_bp, bf$pos_slope_after_bp)
  expect_identical(cs$inside_ci, bf$inside_ci)
})

test_that("engines return NULL below 5 rows", {
  d <- test_data()[1:4, ]
  expect_null(piecewise_metrics_cumsum(d, .x = "vco2", .y = "ve"))
  expect_null(piecewise_metrics_cumsum_jm(d, .x = "vco2", .y = "ve"))
})

test_that("breakpoint() pipeline is stable (v-slope VT1, JM VT2)", {
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
  bp <- breakpoint(
    avg,
    algorithm_vt2 = "jm", x_vt2 = "vco2", y_vt2 = "ve",
    algorithm_vt1 = "v-slope", x_vt1 = "vo2", y_vt1 = "vco2",
    vo2 = "vo2", vco2 = "vco2", ve = "ve", time = "time",
    truncate = TRUE, plots = FALSE
  )

  expect_true(all(bp$bp_dat$determinant_bp))
  vt1 <- bp$bp_dat[bp$bp_dat$bp == "vt1", ]
  vt2 <- bp$bp_dat[bp$bp_dat$bp == "vt2", ]
  expect_equal(unname(vt1$vo2), 2713.058, tolerance = 1e-3)
  expect_equal(unname(vt2$vco2), 2836.667, tolerance = 1e-3)
})
