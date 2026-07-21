# =============================================================================
# Cumulative sufficient statistics — the O(n) piecewise-regression search engine
# =============================================================================
#
# This is the search engine behind loop_orr() and loop_v_slope(): it scores the
# brute-force piecewise-regression search (every division of the data into two
# contiguous OLS segments) in O(n) total time using *cumulative sufficient
# statistics* instead of fitting n separate lm() models. It is mathematically
# identical to the lm()-per-split approach (validated to ~13 significant figures,
# with identical breakpoint selection) but ~250x faster.
#
# It implements the ORR / v-slope "flavor": both regressions are unconstrained
# (ordinary least squares) and the two halves share the split row (left = 1:i,
# right = i:n). It also computes the v-slope-specific distance-to-MSE ratio;
# loop_orr() simply drops that column. The JM flavor (constrained right
# regression) uses a different loop and is not handled here.
#
# The heavy comments below are kept deliberately: this file doubles as the
# readable reference for how the technique works. See
# code_testing/thresholds/piecewise/cumsum_walkthrough.R to step through it.
#
# ---------------------------------------------------------------------------
# THE IDEA
# ---------------------------------------------------------------------------
# For a simple OLS regression y ~ x over any set of points, the slope,
# intercept, and residual sum of squares (RSS) are *exact* functions of just
# six running totals — the "sufficient statistics":
#
#     n,  Sx = sum(x),  Sy = sum(y),  Sxx = sum(x^2),  Syy = sum(y^2),  Sxy = sum(x*y)
#
# From these:
#     Sxx_c = Sxx - Sx^2 / n          # sum of squared deviations of x
#     Sxy_c = Sxy - Sx*Sy / n         # sum of cross-deviations
#     Syy_c = Syy - Sy^2 / n          # sum of squared deviations of y
#     slope     b1 = Sxy_c / Sxx_c
#     intercept b0 = (Sy - b1 * Sx) / n
#     RSS          = Syy_c - Sxy_c^2 / Sxx_c
#
# The trick: to get these totals for EVERY left half (rows 1:i) and EVERY right
# half (rows i:n) at once, precompute prefix sums with cumsum(). Then each
# half's totals are a single lookup / subtraction — O(1) per split, O(n) for
# the whole search. No per-split matrix decomposition.
#
# NUMERICAL NOTE: computing Sxx_c = Sxx - Sx^2/n directly can lose precision
# when x has large magnitude and small spread (exactly like VO2 values in the
# thousands). Subtracting two large near-equal numbers is "catastrophic
# cancellation". Because OLS slope/intercept/RSS are all invariant to shifting
# x and y by a constant, we first CENTER x and y by their global means. That
# shrinks every running total's magnitude and makes the subtraction stable,
# while leaving the slope, RSS, F-test, and intersection unchanged (we shift
# the intersection x back to the original scale at the very end).

#' Piecewise-regression metrics for every split, via cumulative sufficient
#' statistics (O(n)).
#'
#' Shared search engine for `loop_orr()` and `loop_v_slope()`: returns the
#' per-split quantities both need, computed without fitting any lm(). ORR drops
#' the `dist_MSE_ratio` column; v-slope keeps it and renames the direction
#' columns.
#'
#' @param .data A data frame ordered by the x variable.
#' @param .x,.y Column names (character).
#' @param conf_level Confidence level for the approximate JM CI flag.
#'
#' @return A tibble with `idx`, `p`, `pos_change`, `pos_slope_after_bp`,
#'   `int_point_x`, `dist_MSE_ratio`, and `inside_ci` — one row per candidate
#'   split. `NULL` if fewer than 5 rows (a two-line, four-parameter fit is
#'   undefined below that).
#' @keywords internal
#' @noRd
piecewise_metrics_cumsum <- function(.data, .x, .y, conf_level = 0.95) {
  x <- .data[[.x]]
  y <- .data[[.y]]
  n <- length(x)
  if (n < 5) {
    return(NULL)
  }

  # --- center for numerical stability (shift-invariant; undone for x_int) ----
  x_bar <- mean(x)
  y_bar <- mean(y)
  xc <- x - x_bar
  yc <- y - y_bar

  # --- prefix sums with a leading 0 so P[k + 1] = sum of the first k values ---
  # left half 1:i     -> total = P[i + 1]
  # right half i:n    -> total = P[n + 1] - P[i]   (shares row i with the left)
  Px <- c(0, cumsum(xc))
  Py <- c(0, cumsum(yc))
  Pxx <- c(0, cumsum(xc^2))
  Pyy <- c(0, cumsum(yc^2))
  Pxy <- c(0, cumsum(xc * yc))

  i <- seq_len(n)

  # --- left half (rows 1:i), all splits at once -----------------------------
  nL <- i
  SxL <- Px[i + 1]
  SyL <- Py[i + 1]
  SxxL <- Pxx[i + 1]
  SyyL <- Pyy[i + 1]
  SxyL <- Pxy[i + 1]

  # --- right half (rows i:n), all splits at once ----------------------------
  nR <- n - i + 1
  SxR <- Px[n + 1] - Px[i]
  SyR <- Py[n + 1] - Py[i]
  SxxR <- Pxx[n + 1] - Pxx[i]
  SyyR <- Pyy[n + 1] - Pyy[i]
  SxyR <- Pxy[n + 1] - Pxy[i]

  # --- closed-form OLS per half (vectorized over splits) --------------------
  reg <- function(nn, Sx, Sy, Sxx, Syy, Sxy) {
    Sxx_c <- Sxx - Sx^2 / nn
    Sxy_c <- Sxy - Sx * Sy / nn
    Syy_c <- Syy - Sy^2 / nn
    b1 <- Sxy_c / Sxx_c
    b0 <- (Sy - b1 * Sx) / nn
    rss <- Syy_c - Sxy_c^2 / Sxx_c
    list(b0 = b0, b1 = b1, rss = rss)
  }
  L <- reg(nL, SxL, SyL, SxxL, SyyL, SxyL)
  R <- reg(nR, SxR, SyR, SxxR, SyyR, SxyR)

  # --- single-line (simple) model over all n, for the F-test ----------------
  # the n + 1 IS the full model
  Sxx_c_all <- Pxx[n + 1] - Px[n + 1]^2 / n
  Sxy_c_all <- Pxy[n + 1] - Px[n + 1] * Py[n + 1] / n
  Syy_c_all <- Pyy[n + 1] - Py[n + 1]^2 / n
  RSS_simple <- Syy_c_all - Sxy_c_all^2 / Sxx_c_all
  b1_simple <- Sxy_c_all / Sxx_c_all # slope of the single line (centered = raw)

  # --- combine the two halves; F-test on the extra-sum-of-squares -----------
  RSS_two <- L$rss + R$rss
  MSE_two <- RSS_two / (n - 4)
  f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
  p <- stats::pf(f_stat, df1 = 2, df2 = n - 4, lower.tail = FALSE)

  # slope-change direction
  pct_slope_change <- 100 * (R$b1 - L$b1) / abs(L$b1)
  pos_change <- pct_slope_change > 0
  pos_slope_after_bp <- R$b1 > 0

  # intersection of the two lines, in CENTERED coordinates first.
  # matches intersection_point(): x = (b0_right - b0_left) / (b1_left - b1_right)
  int_x_c <- (R$b0 - L$b0) / (L$b1 - R$b1)
  int_y_c <- L$b0 + L$b1 * int_x_c
  int_point_x <- int_x_c + x_bar # shift x back to the original scale

  # v-slope flavor: perpendicular distance from the intersection to the single
  # regression line, over MSE. Distance is translation-invariant, so we compute
  # it in centered coordinates, where the single line passes through the origin
  # (yc = b1_simple * xc). ORR ignores this column; v-slope uses it to rank.
  dist <- abs(b1_simple * int_x_c - int_y_c) / sqrt(b1_simple^2 + 1)
  dist_MSE_ratio <- dist / MSE_two

  # the first and last splits leave a single point in one half -> no line
  boundary <- i == 1 | i == n
  p[boundary] <- NA
  pos_change[boundary] <- NA
  pos_slope_after_bp[boundary] <- NA
  int_point_x[boundary] <- NA
  dist_MSE_ratio[boundary] <- NA
  RSS_two[boundary] <- NA
  MSE_two[boundary] <- NA

  # approximate JM confidence interval: keep splits whose RSS rises less than
  # crit_F * MSE above the minimum RSS.
  crit_F <- stats::qf(conf_level, 1, n - 4, lower.tail = TRUE)
  inside_ci <- dplyr::if_else(
    (RSS_two - min(RSS_two, na.rm = TRUE)) / min(MSE_two, na.rm = TRUE) < crit_F,
    TRUE, FALSE
  )

  tibble::tibble(
    idx = i,
    p = p,
    pos_change = pos_change,
    pos_slope_after_bp = pos_slope_after_bp,
    int_point_x = int_point_x,
    dist_MSE_ratio = dist_MSE_ratio,
    inside_ci = inside_ci
  )
}

piecewise_metrics_cumsum_jm <- function(.data, .x, .y, conf_level = 0.95) {
  x <- .data[[.x]]
  y <- .data[[.y]]
  n <- length(x)
  if (n < 5) { # need more than 5 points for a 4 parameter model
    return(NULL)
  }

  # --- center for numerical stability (shift-invariant; undone for x_int) ----
  x_bar <- mean(x)
  y_bar <- mean(y)
  xc <- x - x_bar
  yc <- y - y_bar

  # --- prefix sums with a leading 0 so P[k + 1] = sum of the first k values ---
  # the leading 0 helps with subsetting the right side
  # left half 1:i     -> total = P[i + 1]
  # JM does not share a split point, so subtract P[i + 1] from P[n + 1]
  # right half i:n    -> total = P[n + 1] - P[i + 1]
  Px <- c(0, cumsum(xc))
  Py <- c(0, cumsum(yc))
  Pxx <- c(0, cumsum(xc^2))
  Pyy <- c(0, cumsum(yc^2))
  Pxy <- c(0, cumsum(xc * yc))

  i <- seq_len(n) # i is the split index

  # --- left half (rows 1:i), all splits at once -----------------------------
  nL <- i # length of left is up to the ith value
  SxL <- Px[i + 1] # add one to account for leading 0
  SyL <- Py[i + 1]
  SxxL <- Pxx[i + 1]
  SyyL <- Pyy[i + 1]
  SxyL <- Pxy[i + 1]

  # using the sufficient statistics for the left side, calculate y0
  # --- closed-form OLS per half (vectorized over splits) --------------------
  reg <- function(nn, Sx, Sy, Sxx, Syy, Sxy) {
    Sxx_c <- Sxx - Sx^2 / nn
    Sxy_c <- Sxy - Sx * Sy / nn
    Syy_c <- Syy - Sy^2 / nn
    b1 <- Sxy_c / Sxx_c
    b0 <- (Sy - b1 * Sx) / nn
    rss <- Syy_c - Sxy_c^2 / Sxx_c
    list(b0 = b0, b1 = b1, rss = rss)
  }
  L <- reg(nL, SxL, SyL, SxxL, SyyL, SxyL)

  # calculate y values at the intersection for all x values (x0 per JM)
  x0 <- xc
  y0 <- L$b0 + L$b1*x0

  # Finding the right half of the JM regression line is akin to shifting the
  # axis to the origin (0,0). Achieve this by subtracting the split point
  # coordinates, x0, and y0, from each x and y value.

  # --- right half (rows i:n), all splits at once ----------------------------
  nR <- n - i # Orr was n - i + 1. The + 1 accounts for the shared split point.

  SxR <- Px[n + 1] - Px[i + 1]
  SyR <- Py[n + 1] - Py[i + 1]
  SxxR <- Pxx[n + 1] - Pxx[i + 1]
  SyyR <- Pyy[n + 1] - Pyy[i + 1]
  SxyR <- Pxy[n + 1] - Pxy[i + 1]

  # x0 per the JM paper is all the x values that you could split at, which is
  # the same thing as xc.
  x0 <- xc

  u = xc - x0
  v = yc - y0

  # Suu = Σ(x − x0)² = SxxR − 2·x0·SxR + nR·x0² (Suu is shifted by x0)
  # Σ(x - x0)(x - x0)
  # Σ(x^2 - 2 * x * x0 + x0^2)
  # Σx^2 - 2 * Σx * x0 + Σx0^2).
  # Σ is equivalent to multiplying by n
  # Σx^2 is SxxR, Σx is SxR
  Suu <- SxxR - 2 * x0 * SxR + nR * x0^2
  Svv <- SyyR - 2 * y0 * SyR + nR * y0^2
  Suv <- SxyR - x0 * SyR - y0 * SxR + nR * x0 * y0

  # because we've centered through the origin, we won't use reg()
  # calculate slope of right regression line
  b3 <- Suv / Suu
  # calculate (centered?) intercept of right regression line
  b2 <- L$b0 + L$b1*x0
  RSS_right <- Svv - Suv^2 / Suu

  # --- single-line (simple) model over all n, for the F-test ----------------
  # the n + 1 IS the full model
  Sxx_c_all <- Pxx[n + 1] - Px[n + 1]^2 / n
  Sxy_c_all <- Pxy[n + 1] - Px[n + 1] * Py[n + 1] / n
  Syy_c_all <- Pyy[n + 1] - Py[n + 1]^2 / n
  RSS_simple <- Syy_c_all - Sxy_c_all^2 / Sxx_c_all
  b1_simple <- Sxy_c_all / Sxx_c_all # slope of the single line (centered = raw)

  # --- combine the two halves; F-test on the extra-sum-of-squares -----------
  RSS_two <- L$rss + RSS_right
  MSE_two <- RSS_two / (n - 4)
  f_stat <- (RSS_simple - RSS_two) / (2 * MSE_two)
  p <- stats::pf(f_stat, df1 = 2, df2 = n - 4, lower.tail = FALSE)

  # slope-change direction
  pct_slope_change <- 100 * (b3 - L$b1) / abs(L$b1)
  pos_change <- pct_slope_change > 0
  pos_slope_after_bp <- b3 > 0

  # intersection of the two lines, in CENTERED coordinates first.
  # matches intersection_point(): x = (b0_right - b0_left) / (b1_left - b1_right)
  int_x_c <- (b2 - L$b0) / (L$b1 - b3)
  int_y_c <- L$b0 + L$b1 * int_x_c
  int_point_x <- int_x_c + x_bar # shift x back to the original scale

  # the first and last splits leave a single point in one half -> no line
  boundary <- i == 1 | i == n
  p[boundary] <- NA
  pos_change[boundary] <- NA
  pos_slope_after_bp[boundary] <- NA
  int_point_x[boundary] <- NA
  dist_MSE_ratio[boundary] <- NA
  RSS_two[boundary] <- NA
  MSE_two[boundary] <- NA

  # approximate JM confidence interval: keep splits whose RSS rises less than
  # crit_F * MSE above the minimum RSS.
  crit_F <- stats::qf(conf_level, 1, n - 4, lower.tail = TRUE)
  inside_ci <- dplyr::if_else(
    (RSS_two - min(RSS_two, na.rm = TRUE)) / min(MSE_two, na.rm = TRUE) < crit_F,
    TRUE, FALSE
  )

  tibble::tibble(
    idx = i,
    p = p,
    pos_change = pos_change,
    pos_slope_after_bp = pos_slope_after_bp,
    int_point_x = int_point_x,
    dist_MSE_ratio = dist_MSE_ratio,
    inside_ci = inside_ci
  )
}
