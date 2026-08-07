# Shared nonparametric bootstrap engine for the derivative-based threshold
# methods (d2_poly_reg_maxima, d2_reg_spline_maxima, d1_crossing, d2_inflection).
#
# Each algorithm supplies its own `statistic` function (case-resampling; see
# boot::boot()'s `statistic` argument) that refits the model on a resampled
# data set and returns a single x-location, or NA_real_ if no valid threshold
# could be found in that resample. This engine handles the parts that are
# identical across algorithms: running boot::boot(), and turning the replicate
# distribution into a percentile confidence interval WITHOUT silently
# producing garbage when some replicates fail (a real bug found in the
# original implementation: on a 27-row resampled data set, 17% of replicates
# returned NA, and broom::tidy(..., conf.method = "perc") did not guard
# against this, yielding an NA point estimate and a nonsensical interval).
#
# Resampling method: ordinary (IID) case resampling of rows, matching the
# rest of this package. NOTE: the fitted curves are typically rolling
# averages of breath-by-breath data and are therefore strongly
# autocorrelated (lag-1 autocorrelation around 0.97 is common). IID
# resampling does not preserve that serial structure, so the resulting
# interval likely UNDERSTATES the true sampling variability of the
# threshold. A block/moving-block bootstrap would be more appropriate for
# this kind of data, but that is a larger methodological change and is not
# implemented here.
#
# No default parallelization: boot::boot(parallel = "multicore") only uses
# more than one core if `ncpus` is also supplied (it defaults to
# getOption("boot.ncpus", 1L) i.e. serial), and "multicore" is unsupported on
# Windows. Since it provided no benefit as configured in the original code,
# it is omitted entirely here rather than carried forward as a no-op.

#' Bootstrap a percentile confidence interval for a threshold location.
#'
#' @param .data Data frame to resample (already front-trimmed / ordered).
#' @param statistic A function with signature `function(.data, i, ...)`
#'   suitable for [boot::boot()]: given the resampled row indices `i`, it must
#'   return a single numeric x-location, or `NA_real_` if no valid threshold
#'   was found in that resample.
#' @param ... Additional arguments passed to `statistic` on every call.
#' @param R Number of bootstrap replicates. Default \code{999}, per standard
#'   guidance for stable percentile confidence intervals (R >= 1000 is
#'   typically recommended; 999 plus the observed statistic gives 1000 values).
#' @param conf_level Confidence level for the interval.
#' @param min_valid Minimum number of non-NA replicates required to report a
#'   confidence interval at all. Below this, the interval is considered too
#'   unstable to report.
#' @param min_valid_frac Minimum fraction of non-NA replicates required (in
#'   addition to `min_valid`).
#'
#' @returns A list with `estimate` (the statistic on the original data),
#'   `lower`, `upper`, `n_valid`, and `n_total`; or `NULL` if too few
#'   replicates produced a valid threshold to estimate an interval (a warning
#'   is issued explaining why).
#' @keywords internal
#' @noRd
bootstrap_threshold_ci <- function(.data, statistic, ...,
                                   R = 999,
                                   conf_level = 0.95,
                                   min_valid = 30,
                                   min_valid_frac = 0.5) {
  boot_res <- boot::boot(
    data = .data, statistic = statistic, R = R, sim = "ordinary", ...
  )

  replicates <- boot_res$t[, 1]
  n_total <- length(replicates)
  n_valid <- sum(!is.na(replicates))

  if (n_valid < min_valid || n_valid / n_total < min_valid_frac) {
    rlang::warn(paste0(
      "Bootstrap confidence interval could not be reliably estimated: only ",
      n_valid, " of ", n_total, " resamples produced a valid threshold. ",
      "Returning the point estimate without a confidence interval."
    ))
    return(NULL)
  }

  if (n_valid < n_total) {
    rlang::warn(paste0(
      n_total - n_valid, " of ", n_total, " bootstrap resamples did not ",
      "produce a valid threshold and were excluded from the confidence ",
      "interval calculation."
    ))
  }

  alpha <- 1 - conf_level
  quantiles <- stats::quantile(
    replicates,
    probs = c(alpha / 2, 1 - alpha / 2),
    na.rm = TRUE
  )

  list(
    estimate = boot_res$t0,
    lower = unname(quantiles[1]),
    upper = unname(quantiles[2]),
    n_valid = n_valid,
    n_total = n_total
  )
}
