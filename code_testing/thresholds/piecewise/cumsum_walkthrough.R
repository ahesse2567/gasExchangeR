# =============================================================================
# Walkthrough: step through piecewise_metrics_cumsum() in the debugger
# =============================================================================
#
# Goal: watch the cumulative-sufficient-statistics regression compute every
# split's slope, intercept, RSS, F-test, and CI flag WITHOUT fitting a single
# lm(), then confirm it matches the brute-force loop.
#
# HOW TO DRIVE THE DEBUGGER (RStudio or console):
#   Run the setup block first (through the `d <- ...` line). Then run the
#   `debugonce(...)` line and the call right after it. R drops you at the top
#   of piecewise_metrics_cumsum(). At the `Browse>` prompt:
#     n        run the next line
#     <Enter>  repeat the last `n` (handy for stepping quickly)
#     s        step INTO a function call on the current line (e.g. reg())
#     f        finish the current function and return
#     c        continue to the end (or next browser)
#     Q        quit debugging
#     where    show the call stack
#   Type ANY variable name at the prompt to print it, e.g. `Px`, `head(Px)`,
#   `L$rss`, `nL[1:5]`. You are in the function's environment, so every local
#   is visible as it comes into existence.
#
# WHAT TO LOOK AT, STAGE BY STAGE (see the "INSPECT" notes below).
# =============================================================================

library(devtools)

# If your working directory is the package root, `.` works. Otherwise point
# load_all() at the gasExchangeR folder.
load_all(".")

# --- 1. Build a realistic dataset through the normal pipeline ---------------
cpet_bbb <- utils::read.csv(
  system.file("extdata", "anton_vo2max_clean.csv", package = "gasExchangeR")
)

cpet_clean <- ventilatory_outliers(
  cpet_bbb,
  outlier_cols = "vo2", time = "time",
  sd_lim = 4, width = 5, mos = "mean", align = "center",
  use_global_sd = TRUE, global_sd_mos = "median",
  exclude_test_val = TRUE, remove_outliers = TRUE,
  max_passes = 1, plot_outliers = FALSE
)

cpet_avg <- avg_exercise_test(
  cpet_clean,
  test = "breath", calc_type = "rolling",
  time_col = "time", window = 15
)

# We'll find VT2 (respiratory compensation): VE on the y-axis, VCO2 on the x.
d <- cpet_avg
x_var <- "vco2"
y_var <- "ve"

# =============================================================================
# 2. Step through the function
# =============================================================================
# debugonce() arms the debugger for exactly one call, without editing source.
debugonce(gasExchangeR:::piecewise_metrics_cumsum)

result <- gasExchangeR:::piecewise_metrics_cumsum(
  d,
  .x = x_var,
  .y = y_var,
  conf_level = 0.95
)

# ---------------------------------------------------------------------------
# INSPECT as you step (type these at the Browse> prompt):
#
#  After `x <- .data[[.x]]` / `y <- ...`:
#     length(x); head(x); n
#
#  After the centering block (xc, yc):
#     mean(xc)            # ~0 by construction
#     range(x); range(xc) # note how much smaller xc is -> numerical stability
#
#  After the prefix sums (Px, Py, Pxx, Pyy, Pxy):
#     head(Px); length(Px)      # length is n+1, with a leading 0
#     Px[1]                     # the leading 0
#     Px[n + 1]                 # equals sum(xc) (~0 after centering)
#     # verify the range-sum trick for, say, rows 10:20:
#     sum(xc[10:20]); Px[21] - Px[10]   # these two should match
#
#  After the left/right totals (SxL, SxR, ...):
#     nL[1:5]; nR[1:5]          # left grows 1,2,3..., right shrinks n, n-1...
#     nL + nR                   # = n + 1 everywhere (they SHARE the split row)
#
#  Step INTO reg() with `s` on the `L <- reg(...)` line to watch one segment's
#  closed-form OLS get built:
#     Sxx_c; Sxy_c; Syy_c       # centered sums (vectors, one per split)
#     b1[100]; b0[100]          # slope/intercept at split 100
#     rss[100]
#  Then `f` to finish reg() and return to the main function.
#
#  After RSS_two / f_stat / p:
#     head(RSS_two); which.min(RSS_two)   # the min-RSS split = best breakpoint
#     head(p)
#
#  After int_point_x:
#     summary(int_point_x)      # note the "+ x_bar" put it back on real scale
#
#  After the boundary NA block:
#     result_preview <- p       # boundaries (1 and n) are NA
#
#  `c` to finish and return the tibble into `result`.
# ---------------------------------------------------------------------------

# 3. Look at the output
print(result, n = 15)

# best breakpoint by the cumsum method (smallest p / lowest RSS in range)
cs_best <- gasExchangeR:::get_best_piecewise_idx(
  result, range(d[[x_var]]),
  alpha_linearity = 0.05, pos_change = TRUE, pos_slope_after_bp = TRUE
)
cat("cumsum best breakpoint index:", cs_best,
    "-> vco2 =", round(d[[x_var]][cs_best], 1), "\n")

# =============================================================================
# 4. Confirm it matches the brute-force lm()-per-split loop
# =============================================================================
bf <- gasExchangeR:::loop_orr(d, .x = x_var, .y = y_var, conf_level = 0.95)
cs <- result[, names(bf)] # same column order

ok <- !is.na(bf$p) # ignore the NA boundary rows
cat("\n--- cumsum vs brute force ---\n")
cat("max |p diff|:          ", format(max(abs(bf$p[ok] - cs$p[ok])), digits = 3), "\n")
cat("max |x_int diff|:      ", format(max(abs(bf$int_point_x[ok] - cs$int_point_x[ok])), digits = 3), "\n")
cat("pos_change identical:  ", identical(bf$pos_change, cs$pos_change), "\n")
cat("inside_ci identical:   ", identical(bf$inside_ci, cs$inside_ci), "\n")

bf_best <- gasExchangeR:::get_best_piecewise_idx(
  bf, range(d[[x_var]]),
  alpha_linearity = 0.05, pos_change = TRUE, pos_slope_after_bp = TRUE
)
cat("same breakpoint chosen:", bf_best == cs_best,
    sprintf("(brute=%d, cumsum=%d)\n", bf_best, cs_best))

# =============================================================================
# 5. (Optional) SEE catastrophic cancellation for yourself
# =============================================================================
# Re-derive one segment's RSS by hand, WITHOUT centering, on x shifted way up.
# Watch the naive formula drift from lm()'s RSS as the magnitudes explode.
demo_cancellation <- function(shift) {
  xx <- d[[x_var]] + shift
  yy <- d[[y_var]]
  n <- length(xx)
  Sx <- sum(xx); Sy <- sum(yy)
  Sxx <- sum(xx^2); Syy <- sum(yy^2); Sxy <- sum(xx * yy)
  Sxx_c <- Sxx - Sx^2 / n
  Sxy_c <- Sxy - Sx * Sy / n
  Syy_c <- Syy - Sy^2 / n
  rss_naive <- Syy_c - Sxy_c^2 / Sxx_c
  rss_lm <- sum(resid(lm(yy ~ xx))^2)
  c(shift = shift, rss_naive = rss_naive, rss_lm = rss_lm,
    abs_diff = abs(rss_naive - rss_lm))
}
cat("\n--- naive (uncentered) RSS vs lm() as x is shifted up ---\n")
print(t(sapply(c(0, 1e3, 1e6, 1e9), demo_cancellation)), digits = 10)
# The abs_diff grows as the shift grows: that is why piecewise_metrics_cumsum()
# centers x and y first.
